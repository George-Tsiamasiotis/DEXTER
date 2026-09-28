//! Custom Error types.

/// Simulation Error type.
#[derive(thiserror::Error, Debug)]
pub enum SimulationError {
    /// From [`dexter_machine::EvalError`].
    #[error("{0}")]
    EvalError(#[from] dexter_machine::EvalError),

    /// From [`dexter_machine::MachineError`].
    #[error("{0}")]
    MachineError(#[from] dexter_machine::MachineError),

    /// Queue initial conditions arrays must have a length of at least 1.
    #[error("Queue initial conditions arrays must have a length of at least 1")]
    QueueInitialConditionsEmptyInput,

    /// Queue initial conditions arrays must be of the same size.
    #[error("Queue initial conditions arrays must be of the same size")]
    QueueInitialConditionsMismatch,
}
