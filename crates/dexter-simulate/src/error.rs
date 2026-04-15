//! Custom Error types.

/// Simulation Error type.
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
}

/// Constants of Motion calculations Error type.
#[derive(thiserror::Error, Debug)]
pub enum COMError {
    /// From [`dexter_equilibrium::EvalError`].
    #[error("{0}")]
    EvalError(#[from] dexter_equilibrium::EvalError),

    /// From [`dexter_equilibrium::EqError`].
    #[error("{0}")]
    EqError(#[from] dexter_equilibrium::EqError),

    /// [`COMs`](crate::COMs) missing 'mu' field.
    #[error("'COMs' missing 'mu' field")]
    UndefinedMu,

    /// [`COMs`](crate::COMs) missing 'pzeta' field.
    #[error("'COMs' missing 'pzeta' field")]
    UndefinedPzeta,

    /// [`COMs`](crate::COMs) missing 'energy' field.
    #[error("'COMs' missing 'energy' field")]
    UndefinedEnergy,
}
