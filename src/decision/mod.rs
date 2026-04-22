pub mod engine;
pub mod negative_control;
pub mod release_gate;
pub mod stats;

pub use engine::{DecisionContext, DecisionEngine, DecisionOutcome};
pub use negative_control::{
    apply_negative_control, BackgroundComparator, NegativeControlDecisionInput,
    NegativeControlResult,
};
pub use release_gate::{load_benchmark_gates, ReleaseContext, ReleaseGate};
pub use stats::{ProportionEstimate, SAMPLING_ONLY_CI_LABEL};
