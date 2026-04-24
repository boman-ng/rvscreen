use super::spool::{RepresentativeBatchDelta, RepresentativeRoundSpool};
use crate::adjudicate::{AdjudicationInputs, FragmentAdjudicator};
use crate::aggregate::CandidateAggregator;
use crate::align::{CompetitiveAligner, FragmentAlignResult};
use crate::error::{Result, RvScreenError};
use crate::io::FragmentRecord as IoFragmentRecord;
use crate::types::FragmentClass;
use std::collections::{BTreeMap, VecDeque};
use std::marker::PhantomData;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc::{sync_channel, Receiver, SyncSender},
    Arc, Condvar, Mutex,
};
use std::thread::{self, ScopedJoinHandle};

const FRAGMENTS_PER_THREAD: usize = 64;
const MIN_ALIGNMENT_BATCH_SIZE: usize = 64;
const MAX_ALIGNMENT_BATCH_SIZE: usize = 1024;
const MAX_EXECUTOR_THREADS: usize = 16;

pub(crate) struct RoundParallelism {
    threads: usize,
    batch_size: usize,
}

struct ProcessedFragment {
    class: FragmentClass,
    align_result: FragmentAlignResult,
}

struct RepresentativeSubmission {
    entry_round: usize,
    fragment: IoFragmentRecord,
}

struct RepresentativeBatch {
    batch_sequence: usize,
    submissions: Vec<RepresentativeSubmission>,
}

struct Sequenced<T> {
    sequence: usize,
    payload: T,
}

struct SharedFailure {
    cancelled: AtomicBool,
    first_message: Mutex<Option<String>>,
}

struct WorkQueue<T> {
    state: Mutex<WorkQueueState<T>>,
    submit_ready: Condvar,
    worker_ready: Vec<Condvar>,
}

struct WorkQueueState<T> {
    queues: Vec<VecDeque<Sequenced<T>>>,
    capacities: Vec<usize>,
    next_worker: usize,
    closed: bool,
    cancelled: bool,
}

struct OrderedExecutor<'scope, In, Out, SinkState>
where
    In: Send + 'scope,
    Out: Send + 'scope,
    SinkState: Send + 'scope,
{
    work_queue: Arc<WorkQueue<In>>,
    worker_handles: Vec<ScopedJoinHandle<'scope, Result<()>>>,
    sink_handle: Option<ScopedJoinHandle<'scope, Result<SinkState>>>,
    failure: Arc<SharedFailure>,
    _out: PhantomData<fn() -> Out>,
}

pub(crate) struct RoundExecutor<'scope> {
    inner: OrderedExecutor<'scope, IoFragmentRecord, ProcessedFragment, CandidateAggregator>,
}

pub(crate) struct RepresentativeRoundExecutor<'scope> {
    inner: OrderedExecutor<
        'scope,
        RepresentativeBatch,
        RepresentativeBatchDelta,
        RepresentativeRoundSpool,
    >,
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
            threads: requested_threads.min(MAX_EXECUTOR_THREADS),
            batch_size: batch_size_for(requested_threads.min(MAX_EXECUTOR_THREADS)),
        })
    }

    pub(crate) fn threads(&self) -> usize {
        self.threads
    }

    pub(crate) fn representative_batch_size(&self) -> usize {
        self.batch_size
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

    pub(crate) fn start_representative_round_executor<'scope>(
        &self,
        scope: &'scope thread::Scope<'scope, '_>,
        aligner: &'scope CompetitiveAligner,
        adjudicator: &'scope FragmentAdjudicator,
        spool: RepresentativeRoundSpool,
    ) -> RepresentativeRoundExecutor<'scope> {
        let round_count = spool.round_count();
        let work_capacity = self.threads.max(1);
        let result_capacity = self.threads.max(1);
        let inner = OrderedExecutor::new(
            scope,
            self.threads,
            work_capacity,
            result_capacity,
            move |batch| process_representative_batch(batch, round_count, aligner, adjudicator),
            spool,
            |spool, batch_delta| spool.merge_batch_delta(batch_delta),
        );

        RepresentativeRoundExecutor { inner }
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

impl<'scope> RepresentativeRoundExecutor<'scope> {
    pub(crate) fn submit_batch(
        &self,
        batch_sequence: usize,
        fragments: Vec<(usize, IoFragmentRecord)>,
    ) -> Result<()> {
        self.inner.submit(
            batch_sequence,
            RepresentativeBatch {
                batch_sequence,
                submissions: fragments
                    .into_iter()
                    .map(|(entry_round, fragment)| RepresentativeSubmission {
                        entry_round,
                        fragment,
                    })
                    .collect(),
            },
        )
    }

    pub(crate) fn finish(self) -> Result<RepresentativeRoundSpool> {
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

impl<T> WorkQueue<T> {
    fn new(worker_count: usize, total_capacity: usize) -> Self {
        let capacities = distribute_work_capacity(worker_count.max(1), total_capacity.max(1));
        let queues = capacities
            .iter()
            .map(|capacity| VecDeque::with_capacity(*capacity))
            .collect();
        let worker_ready = (0..worker_count.max(1)).map(|_| Condvar::new()).collect();

        Self {
            state: Mutex::new(WorkQueueState {
                queues,
                capacities,
                next_worker: 0,
                closed: false,
                cancelled: false,
            }),
            submit_ready: Condvar::new(),
            worker_ready,
        }
    }

    fn submit(&self, item: Sequenced<T>, failure: &SharedFailure) -> Result<()> {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        loop {
            if state.cancelled || failure.is_cancelled() {
                return Err(failure.stage_error(
                    "bounded topology submission aborted",
                    "downstream stage reported a failure while the work queue was full",
                ));
            }

            if state.closed {
                return Err(failure.stage_error(
                    "bounded topology channel closed",
                    "work queue disconnected before submission completed",
                ));
            }

            if let Some(worker_index) = next_available_worker(&state) {
                state.queues[worker_index].push_back(item);
                state.next_worker = (worker_index + 1) % state.queues.len();
                drop(state);
                self.worker_ready[worker_index].notify_one();
                return Ok(());
            }

            state = self
                .submit_ready
                .wait(state)
                .unwrap_or_else(|poisoned| poisoned.into_inner());
        }
    }

    fn recv(&self, worker_index: usize) -> Option<Sequenced<T>> {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());

        loop {
            if state.cancelled {
                return None;
            }

            if let Some(item) = state.queues[worker_index].pop_front() {
                drop(state);
                self.submit_ready.notify_one();
                return Some(item);
            }

            if state.closed {
                return None;
            }

            state = self.worker_ready[worker_index]
                .wait(state)
                .unwrap_or_else(|poisoned| poisoned.into_inner());
        }
    }

    fn close(&self) {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if state.closed {
            return;
        }
        state.closed = true;
        drop(state);
        self.submit_ready.notify_all();
        self.notify_workers();
    }

    fn cancel(&self) {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if state.cancelled {
            return;
        }
        state.cancelled = true;
        for queue in &mut state.queues {
            queue.clear();
        }
        drop(state);
        self.submit_ready.notify_all();
        self.notify_workers();
    }

    fn notify_workers(&self) {
        for worker_ready in &self.worker_ready {
            worker_ready.notify_all();
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
        let worker_count = threads.clamp(1, MAX_EXECUTOR_THREADS);
        let work_queue = Arc::new(WorkQueue::new(worker_count, work_capacity));
        let (result_tx, result_rx) = sync_channel(result_capacity.max(1));
        let failure = Arc::new(SharedFailure::new());
        let worker_fn = Arc::new(worker_fn);
        let aggregate_fn = Arc::new(aggregate_fn);

        let mut worker_handles = Vec::with_capacity(worker_count);
        for worker_index in 0..worker_count {
            let work_queue = Arc::clone(&work_queue);
            let result_tx = result_tx.clone();
            let worker_fn = Arc::clone(&worker_fn);
            let failure = Arc::clone(&failure);
            worker_handles.push(scope.spawn(move || {
                worker_loop(worker_index, work_queue, result_tx, worker_fn, failure)
            }));
        }
        drop(result_tx);

        let work_queue_for_sink = Arc::clone(&work_queue);
        let failure_for_sink = Arc::clone(&failure);
        let sink_handle = scope.spawn(move || {
            sink_loop(
                result_rx,
                sink_state,
                aggregate_fn,
                work_queue_for_sink,
                failure_for_sink,
            )
        });

        Self {
            work_queue,
            worker_handles,
            sink_handle: Some(sink_handle),
            failure,
            _out: PhantomData,
        }
    }

    fn submit(&self, sequence: usize, payload: In) -> Result<()> {
        self.work_queue
            .submit(Sequenced { sequence, payload }, self.failure.as_ref())
    }

    fn finish(mut self) -> Result<SinkState> {
        self.work_queue.close();

        let mut first_error = None;
        for handle in std::mem::take(&mut self.worker_handles) {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) => {
                    self.work_queue.cancel();
                    if first_error.is_none() {
                        first_error = Some(error);
                    }
                }
                Err(_) => {
                    self.work_queue.cancel();
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

impl<'scope, In, Out, SinkState> Drop for OrderedExecutor<'scope, In, Out, SinkState>
where
    In: Send + 'scope,
    Out: Send + 'scope,
    SinkState: Send + 'scope,
{
    fn drop(&mut self) {
        self.work_queue.cancel();
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

fn process_representative_batch(
    batch: RepresentativeBatch,
    round_count: usize,
    aligner: &CompetitiveAligner,
    adjudicator: &FragmentAdjudicator,
) -> Result<RepresentativeBatchDelta> {
    let mut batch_delta = RepresentativeBatchDelta::new(batch.batch_sequence, round_count);

    for submission in batch.submissions {
        let processed = process_fragment(submission.fragment, aligner, adjudicator)?;
        batch_delta.record_processed_fragment(
            submission.entry_round,
            processed.class,
            &processed.align_result,
        )?;
    }

    Ok(batch_delta)
}

fn batch_size_for(threads: usize) -> usize {
    threads
        .saturating_mul(FRAGMENTS_PER_THREAD)
        .clamp(MIN_ALIGNMENT_BATCH_SIZE, MAX_ALIGNMENT_BATCH_SIZE)
}

fn distribute_work_capacity(worker_count: usize, total_capacity: usize) -> Vec<usize> {
    let worker_count = worker_count.max(1);
    let total_capacity = total_capacity.max(1);
    let base_capacity = total_capacity / worker_count;
    let remainder = total_capacity % worker_count;

    (0..worker_count)
        .map(|worker_index| base_capacity + usize::from(worker_index < remainder))
        .collect()
}

fn next_available_worker<T>(state: &WorkQueueState<T>) -> Option<usize> {
    for offset in 0..state.queues.len() {
        let worker_index = (state.next_worker + offset) % state.queues.len();
        if state.queues[worker_index].len() < state.capacities[worker_index] {
            return Some(worker_index);
        }
    }

    None
}

fn worker_loop<In, Out, WorkerFn>(
    worker_index: usize,
    work_queue: Arc<WorkQueue<In>>,
    result_tx: SyncSender<Sequenced<Out>>,
    worker_fn: Arc<WorkerFn>,
    failure: Arc<SharedFailure>,
) -> Result<()>
where
    In: Send,
    Out: Send,
    WorkerFn: Fn(In) -> Result<Out> + Send + Sync,
{
    while let Some(work_item) = work_queue.recv(worker_index) {
        let output = match worker_fn(work_item.payload) {
            Ok(output) => output,
            Err(error) => {
                failure.record(&error);
                work_queue.cancel();
                return Err(error);
            }
        };

        let result_item = Sequenced {
            sequence: work_item.sequence,
            payload: output,
        };
        if let Err(error) = send_result(
            &result_tx,
            result_item,
            "result queue disconnected before worker result could be emitted",
            failure.as_ref(),
        ) {
            failure.record(&error);
            work_queue.cancel();
            return Err(error);
        }
    }

    Ok(())
}

fn sink_loop<In, Out, SinkState, AggregateFn>(
    result_rx: Receiver<Sequenced<Out>>,
    mut sink_state: SinkState,
    aggregate_fn: Arc<AggregateFn>,
    work_queue: Arc<WorkQueue<In>>,
    failure: Arc<SharedFailure>,
) -> Result<SinkState>
where
    In: Send,
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
            work_queue.cancel();
            return Err(error);
        }
        if let Err(error) = drain_pending_results(
            &mut pending,
            &mut next_sequence,
            &mut sink_state,
            aggregate_fn.as_ref(),
            failure.as_ref(),
        ) {
            work_queue.cancel();
            return Err(error);
        }
    }

    if failure.is_cancelled() {
        work_queue.cancel();
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
        work_queue.cancel();
        return Err(error);
    }

    Ok(sink_state)
}

fn drain_pending_results<Out, SinkState, AggregateFn>(
    pending: &mut BTreeMap<usize, Out>,
    next_sequence: &mut usize,
    sink_state: &mut SinkState,
    aggregate_fn: &AggregateFn,
    failure: &SharedFailure,
) -> Result<()>
where
    AggregateFn: Fn(&mut SinkState, Out) -> Result<()>,
{
    while let Some(item) = pending.remove(next_sequence) {
        if let Err(error) = aggregate_fn(sink_state, item) {
            failure.record(&error);
            return Err(error);
        }
        *next_sequence = next_sequence.saturating_add(1);
    }

    Ok(())
}

fn send_result<T>(
    sender: &SyncSender<Sequenced<T>>,
    item: Sequenced<T>,
    disconnected_reason: &str,
    failure: &SharedFailure,
) -> Result<()> {
    sender
        .send(item)
        .map_err(|_| failure.stage_error("bounded topology channel closed", disconnected_reason))
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

    #[test]
    fn ordered_executor_preserves_submission_order_across_worker_reordering() {
        thread::scope(|scope| {
            let executor = OrderedExecutor::new(
                scope,
                4,
                4,
                1,
                |value: usize| {
                    let delay_ms = ((8usize.saturating_sub(value)) * 10) as u64;
                    std::thread::sleep(Duration::from_millis(delay_ms));
                    Ok(value)
                },
                Vec::new(),
                |seen: &mut Vec<usize>, value| {
                    seen.push(value);
                    Ok(())
                },
            );

            for value in 0..8 {
                executor
                    .submit(value, value)
                    .expect("executor should accept ordered work");
            }

            let seen = executor.finish().expect("executor should finish cleanly");
            assert_eq!(seen, (0..8).collect::<Vec<_>>());
        });
    }

    #[test]
    fn ordered_executor_unblocks_submitter_after_worker_failure() {
        let (started_tx, started_rx) = sync_channel(1);
        let (release_tx, release_rx) = sync_channel(1);
        let release_rx = Arc::new(Mutex::new(release_rx));

        thread::scope(|scope| {
            let release_rx = Arc::clone(&release_rx);
            let executor = OrderedExecutor::new(
                scope,
                1,
                1,
                1,
                move |value: usize| {
                    if value == 0 {
                        started_tx
                            .send(())
                            .expect("worker should signal when it starts processing");
                        let receiver = release_rx
                            .lock()
                            .unwrap_or_else(|poisoned| poisoned.into_inner());
                        receiver
                            .recv()
                            .expect("worker release signal should arrive before failure");
                        return Err(RvScreenError::validation(
                            "pipeline.parallel",
                            "intentional worker failure",
                        ));
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
                .expect("second item should occupy the only queued slot");
            started_rx
                .recv_timeout(Duration::from_secs(1))
                .expect("worker should start processing the failing item");

            let release_handle = std::thread::spawn(move || {
                std::thread::sleep(Duration::from_millis(150));
                release_tx
                    .send(())
                    .expect("release signal should unblock the failing worker");
            });

            let started_wait = Instant::now();
            let submit_error = executor
                .submit(2, 2)
                .expect_err("blocked submit should return once worker failure cancels the queue");
            assert!(
                started_wait.elapsed() >= Duration::from_millis(100),
                "blocked submit should wait until the worker failure is observed"
            );
            assert!(
                submit_error
                    .to_string()
                    .contains("intentional worker failure"),
                "submit error should include the original worker failure: {submit_error}"
            );

            let finish_error = executor
                .finish()
                .expect_err("finish should surface the worker failure");
            assert!(
                finish_error
                    .to_string()
                    .contains("intentional worker failure"),
                "finish error should include the original worker failure: {finish_error}"
            );

            release_handle
                .join()
                .expect("release helper thread should not panic");
        });
    }

    #[test]
    fn round_parallelism_clamps_executor_threads_to_thread_budget() {
        let parallelism = RoundParallelism::new(MAX_EXECUTOR_THREADS + 8)
            .expect("parallelism should clamp oversized thread requests");

        assert_eq!(parallelism.threads, MAX_EXECUTOR_THREADS);
    }
}
