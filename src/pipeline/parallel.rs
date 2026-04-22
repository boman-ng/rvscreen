use crate::adjudicate::{AdjudicationInputs, FragmentAdjudicator};
use crate::aggregate::CandidateAggregator;
use crate::align::{CompetitiveAligner, FragmentAlignResult};
use crate::error::{Result, RvScreenError};
use crate::io::FragmentRecord as IoFragmentRecord;
use crate::types::FragmentClass;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc::{Receiver, SyncSender, TryRecvError, TrySendError, sync_channel},
    Arc, Mutex,
};
use std::thread::{self, ScopedJoinHandle};

const FRAGMENTS_PER_THREAD: usize = 64;
const MIN_ALIGNMENT_BATCH_SIZE: usize = 64;
const MAX_ALIGNMENT_BATCH_SIZE: usize = 1024;

pub(crate) struct RoundParallelism {
    threads: usize,
    batch_size: usize,
}

struct ProcessedFragment {
    class: FragmentClass,
    align_result: FragmentAlignResult,
}

struct Sequenced<T> {
    sequence: usize,
    payload: T,
}

struct SharedFailure {
    cancelled: AtomicBool,
    first_message: Mutex<Option<String>>,
}

struct OrderedExecutor<'scope, In, Out, SinkState>
where
    In: Send + 'scope,
    Out: Send + 'scope,
    SinkState: Send + 'scope,
{
    work_tx: Option<SyncSender<Sequenced<In>>>,
    submit_permit_tx: SyncSender<()>,
    submit_permit_rx: Receiver<()>,
    worker_handles: Vec<ScopedJoinHandle<'scope, Result<()>>>,
    sink_handle: Option<ScopedJoinHandle<'scope, Result<SinkState>>>,
    failure: Arc<SharedFailure>,
    _out: PhantomData<fn() -> Out>,
}

pub(crate) struct RoundExecutor<'scope> {
    inner: OrderedExecutor<'scope, IoFragmentRecord, ProcessedFragment, CandidateAggregator>,
}

impl RoundParallelism {
    pub(crate) fn new(requested_threads: usize) -> Result<Self> {
        if requested_threads == 0 {
            return Err(RvScreenError::validation(
                "threads",
                "--threads must be greater than zero",
            ));
        }

        Ok(Self {
            threads: requested_threads,
            batch_size: batch_size_for(requested_threads),
        })
    }

    pub(crate) fn start_round_executor<'scope>(
        &self,
        scope: &'scope thread::Scope<'scope, '_>,
        aligner: &'scope CompetitiveAligner,
        adjudicator: &'scope FragmentAdjudicator,
    ) -> RoundExecutor<'scope> {
        let work_capacity = self.batch_size;
        let result_capacity = self.batch_size;
        let inner = OrderedExecutor::new(
            scope,
            self.threads,
            work_capacity,
            result_capacity,
            move |fragment| process_fragment(fragment, aligner, adjudicator),
            CandidateAggregator::new(),
            |aggregator, processed| {
                aggregator.add_fragment(processed.class, &processed.align_result)
            },
        );

        RoundExecutor { inner }
    }
}

impl<'scope> RoundExecutor<'scope> {
    pub(crate) fn submit(&self, sequence: usize, fragment: IoFragmentRecord) -> Result<()> {
        self.inner.submit(sequence, fragment)
    }

    pub(crate) fn finish(self) -> Result<CandidateAggregator> {
        self.inner.finish()
    }
}

impl SharedFailure {
    fn new() -> Self {
        Self {
            cancelled: AtomicBool::new(false),
            first_message: Mutex::new(None),
        }
    }

    fn record(&self, error: &RvScreenError) {
        self.cancelled.store(true, Ordering::Release);
        let mut slot = self
            .first_message
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if slot.is_none() {
            *slot = Some(error.to_string());
        }
    }

    fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Acquire)
    }

    fn stage_error(&self, stage: &str, reason: impl Into<String>) -> RvScreenError {
        let reason = reason.into();
        let message = self
            .first_message
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .clone();
        match message {
            Some(message) => RvScreenError::validation(
                "pipeline.parallel",
                format!("{stage}: {reason}; first failure: {message}"),
            ),
            None => RvScreenError::validation("pipeline.parallel", format!("{stage}: {reason}")),
        }
    }
}

impl<'scope, In, Out, SinkState> OrderedExecutor<'scope, In, Out, SinkState>
where
    In: Send + 'scope,
    Out: Send + 'scope,
    SinkState: Send + 'scope,
{
    fn new<WorkerFn, AggregateFn>(
        scope: &'scope thread::Scope<'scope, '_>,
        threads: usize,
        work_capacity: usize,
        result_capacity: usize,
        worker_fn: WorkerFn,
        sink_state: SinkState,
        aggregate_fn: AggregateFn,
    ) -> Self
    where
        WorkerFn: Fn(In) -> Result<Out> + Send + Sync + 'scope,
        AggregateFn: Fn(&mut SinkState, Out) -> Result<()> + Send + Sync + 'scope,
    {
        let (work_tx, work_rx) = sync_channel(work_capacity.max(1));
        let (result_tx, result_rx) = sync_channel(result_capacity.max(1));
        let max_in_flight = work_capacity.max(1) + result_capacity.max(1) + threads.max(1);
        let (submit_permit_tx, submit_permit_rx) = sync_channel(max_in_flight);
        for _ in 0..max_in_flight {
            submit_permit_tx
                .send(())
                .expect("fresh submit-permit channel should accept its initial capacity");
        }
        let failure = Arc::new(SharedFailure::new());
        let shared_work_rx = Arc::new(Mutex::new(work_rx));
        let worker_fn = Arc::new(worker_fn);
        let aggregate_fn = Arc::new(aggregate_fn);

        let mut worker_handles = Vec::with_capacity(threads.max(1));
        for _ in 0..threads.max(1) {
            let work_rx = Arc::clone(&shared_work_rx);
            let result_tx = result_tx.clone();
            let worker_fn = Arc::clone(&worker_fn);
            let failure = Arc::clone(&failure);
            worker_handles.push(scope.spawn(move || {
                worker_loop(work_rx, result_tx, worker_fn, failure)
            }));
        }
        drop(result_tx);

        let submit_permit_for_sink = submit_permit_tx.clone();
        let failure_for_sink = Arc::clone(&failure);
        let sink_handle = scope.spawn(move || {
            sink_loop(
                result_rx,
                sink_state,
                aggregate_fn,
                submit_permit_for_sink,
                failure_for_sink,
            )
        });

        Self {
            work_tx: Some(work_tx),
            submit_permit_tx,
            submit_permit_rx,
            worker_handles,
            sink_handle: Some(sink_handle),
            failure,
            _out: PhantomData,
        }
    }

    fn submit(&self, sequence: usize, payload: In) -> Result<()> {
        acquire_submit_permit(&self.submit_permit_rx, self.failure.as_ref())?;
        let work_tx = self
            .work_tx
            .as_ref()
            .expect("ordered executor submit called after finish");
        if let Err(error) = send_with_backpressure(
            work_tx,
            Sequenced { sequence, payload },
            "work queue disconnected before submission completed",
            &self.failure,
        ) {
            restore_submit_permit(&self.submit_permit_tx, self.failure.as_ref())?;
            return Err(error);
        }

        Ok(())
    }

    fn finish(mut self) -> Result<SinkState> {
        drop(self.work_tx.take());

        let mut first_error = None;
        for handle in self.worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) => {
                    if first_error.is_none() {
                        first_error = Some(error);
                    }
                }
                Err(_) => {
                    if first_error.is_none() {
                        first_error = Some(RvScreenError::validation(
                            "pipeline.parallel",
                            "round worker panicked",
                        ));
                    }
                }
            }
        }

        let sink_result = match self
            .sink_handle
            .take()
            .expect("ordered executor sink handle should exist")
            .join()
        {
            Ok(Ok(state)) => Ok(state),
            Ok(Err(error)) => Err(error),
            Err(_) => Err(RvScreenError::validation(
                "pipeline.parallel",
                "round aggregation sink panicked",
            )),
        };

        match (first_error, sink_result) {
            (Some(error), _) => Err(error),
            (None, Err(error)) => Err(error),
            (None, Ok(state)) => Ok(state),
        }
    }
}

fn process_fragment(
    fragment: IoFragmentRecord,
    aligner: &CompetitiveAligner,
    adjudicator: &FragmentAdjudicator,
) -> Result<ProcessedFragment> {
    let align_result = aligner.align_fragment(&fragment)?;
    let adjudication = adjudicator.adjudicate(&align_result, AdjudicationInputs::default());

    Ok(ProcessedFragment {
        class: adjudication.class,
        align_result,
    })
}

fn batch_size_for(threads: usize) -> usize {
    threads
        .saturating_mul(FRAGMENTS_PER_THREAD)
        .clamp(MIN_ALIGNMENT_BATCH_SIZE, MAX_ALIGNMENT_BATCH_SIZE)
}

fn worker_loop<In, Out, WorkerFn>(
    work_rx: Arc<Mutex<Receiver<Sequenced<In>>>>,
    result_tx: SyncSender<Sequenced<Out>>,
    worker_fn: Arc<WorkerFn>,
    failure: Arc<SharedFailure>,
) -> Result<()>
where
    In: Send,
    Out: Send,
    WorkerFn: Fn(In) -> Result<Out> + Send + Sync,
{
    loop {
        let work_item = {
            let receiver = work_rx
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            receiver.recv()
        };

        let work_item = match work_item {
            Ok(work_item) => work_item,
            Err(_) => return Ok(()),
        };

        let output = match worker_fn(work_item.payload) {
            Ok(output) => output,
            Err(error) => {
                failure.record(&error);
                return Err(error);
            }
        };

        let result_item = Sequenced {
            sequence: work_item.sequence,
            payload: output,
        };
        if let Err(error) = send_with_backpressure(
            &result_tx,
            result_item,
            "result queue disconnected before worker result could be emitted",
            &failure,
        ) {
            failure.record(&error);
            return Err(error);
        }
    }
}

fn sink_loop<Out, SinkState, AggregateFn>(
    result_rx: Receiver<Sequenced<Out>>,
    mut sink_state: SinkState,
    aggregate_fn: Arc<AggregateFn>,
    submit_permit_tx: SyncSender<()>,
    failure: Arc<SharedFailure>,
) -> Result<SinkState>
where
    Out: Send,
    SinkState: Send,
    AggregateFn: Fn(&mut SinkState, Out) -> Result<()> + Send + Sync,
{
    let mut next_sequence = 0usize;
    let mut pending = BTreeMap::new();

    while let Ok(item) = result_rx.recv() {
        if pending.insert(item.sequence, item.payload).is_some() {
            let error = RvScreenError::validation(
                "pipeline.parallel",
                format!(
                    "received duplicate worker result for sequence {}",
                    item.sequence
                ),
            );
            failure.record(&error);
            return Err(error);
        }
        drain_pending_results(
            &mut pending,
            &mut next_sequence,
            &mut sink_state,
            aggregate_fn.as_ref(),
            &submit_permit_tx,
            failure.as_ref(),
        )?;
    }

    if failure.is_cancelled() {
        return Err(failure.stage_error(
            "aggregation sink terminated after downstream cancellation",
            "round pipeline did not complete cleanly",
        ));
    }

    if !pending.is_empty() {
        let error = failure.stage_error(
            "aggregation sink observed a gap in ordered results",
            format!(
                "missing sequence {next_sequence} while {} later result(s) remained buffered",
                pending.len()
            ),
        );
        failure.record(&error);
        return Err(error);
    }

    Ok(sink_state)
}

fn drain_pending_results<Out, SinkState, AggregateFn>(
    pending: &mut BTreeMap<usize, Out>,
    next_sequence: &mut usize,
    sink_state: &mut SinkState,
    aggregate_fn: &AggregateFn,
    submit_permit_tx: &SyncSender<()>,
    failure: &SharedFailure,
) -> Result<()>
where
    AggregateFn: Fn(&mut SinkState, Out) -> Result<()>,
{
    while let Some(item) = pending.remove(&*next_sequence) {
        if let Err(error) = aggregate_fn(sink_state, item) {
            failure.record(&error);
            return Err(error);
        }
        restore_submit_permit(submit_permit_tx, failure)?;
        *next_sequence = next_sequence.saturating_add(1);
    }

    Ok(())
}

fn acquire_submit_permit(receiver: &Receiver<()>, failure: &SharedFailure) -> Result<()> {
    loop {
        match receiver.try_recv() {
            Ok(()) => return Ok(()),
            Err(TryRecvError::Empty) => {
                if failure.is_cancelled() {
                    return Err(failure.stage_error(
                        "bounded topology submit window closed",
                        "downstream stage reported a failure while waiting for aggregation capacity",
                    ));
                }
                thread::yield_now();
            }
            Err(TryRecvError::Disconnected) => {
                return Err(failure.stage_error(
                    "submit permit channel closed",
                    "aggregation sink stopped issuing capacity permits",
                ));
            }
        }
    }
}

fn restore_submit_permit(sender: &SyncSender<()>, failure: &SharedFailure) -> Result<()> {
    sender.send(()).map_err(|_| {
        failure.stage_error(
            "submit permit release failed",
            "aggregation sink dropped the capacity channel unexpectedly",
        )
    })
}

fn send_with_backpressure<T>(
    sender: &SyncSender<Sequenced<T>>,
    mut item: Sequenced<T>,
    disconnected_reason: &str,
    failure: &Arc<SharedFailure>,
) -> Result<()> {
    loop {
        match sender.try_send(item) {
            Ok(()) => return Ok(()),
            Err(TrySendError::Full(returned)) => {
                if failure.is_cancelled() {
                    return Err(failure.stage_error(
                        "bounded topology submission aborted",
                        "downstream stage reported a failure while the queue was full",
                    ));
                }
                item = returned;
                thread::yield_now();
            }
            Err(TrySendError::Disconnected(_)) => {
                return Err(failure.stage_error("bounded topology channel closed", disconnected_reason));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        atomic::{AtomicBool, Ordering},
        mpsc::sync_channel,
    };
    use std::time::{Duration, Instant};

    #[test]
    fn ordered_executor_applies_backpressure_to_a_full_work_queue() {
        let first_item_gate = Arc::new(AtomicBool::new(true));
        let (started_tx, started_rx) = sync_channel(1);
        let (release_tx, release_rx) = sync_channel(1);
        let release_rx = Arc::new(Mutex::new(release_rx));

        thread::scope(|scope| {
            let worker_gate = Arc::clone(&first_item_gate);
            let release_rx = Arc::clone(&release_rx);
            let executor = OrderedExecutor::new(
                scope,
                1,
                1,
                1,
                move |value: usize| {
                    if worker_gate.swap(false, Ordering::AcqRel) {
                        started_tx
                            .send(())
                            .expect("worker should signal when it starts processing");
                        let receiver = release_rx
                            .lock()
                            .unwrap_or_else(|poisoned| poisoned.into_inner());
                        receiver
                            .recv()
                            .expect("worker gate should be released after backpressure assertion");
                    }
                    Ok(value)
                },
                Vec::new(),
                |seen: &mut Vec<usize>, value| {
                    seen.push(value);
                    Ok(())
                },
            );

            executor
                .submit(0, 0)
                .expect("first item should enter the topology");
            executor
                .submit(1, 1)
                .expect("second item should fill the bounded work queue");
            started_rx
                .recv_timeout(Duration::from_secs(1))
                .expect("worker should begin processing the first item");

            let release_handle = std::thread::spawn(move || {
                std::thread::sleep(Duration::from_millis(150));
                release_tx
                    .send(())
                    .expect("release signal should unblock the worker");
            });

            let started_wait = Instant::now();
            executor
                .submit(2, 2)
                .expect("third item should eventually enqueue after backpressure releases");
            assert!(
                started_wait.elapsed() >= Duration::from_millis(100),
                "third submit should block behind the bounded work queue"
            );

            let seen = executor.finish().expect("executor should finish cleanly");
            release_handle
                .join()
                .expect("release helper thread should not panic");
            assert_eq!(seen, vec![0, 1, 2]);
        });
    }
}
