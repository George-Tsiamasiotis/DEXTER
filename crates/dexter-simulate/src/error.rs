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
}
