//! Custom Error types.

/// Top level Error type.
#[derive(thiserror::Error, Debug)]
pub enum SimulationError {
    /// From [`dexter_equilibrium::EvalError`].
    #[error("{0}")]
    EvalError(#[from] dexter_equilibrium::EvalError),

    /// From [`dexter_equilibrium::EqError`].
    #[error("{0}")]
    EqError(#[from] dexter_equilibrium::EqError),

    /// Queue initial conditions arrays must have a length of at least 1.
    #[error("Queue initial conditions arrays must have a length of at least 1")]
    QueueInitialConditionsEmptyInput,

    /// Queue initial conditions arrays must be of the same size.
    #[error("Queue initial conditions arrays must be of the same size")]
    QueueInitialConditionsMismatch,

    /// Missing required field in [`COMs`](crate::COMs).
    #[error("Missing '{0}' field")]
    MissingCOM(Box<str>),
}
