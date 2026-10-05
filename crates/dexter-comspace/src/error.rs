//! Custom Error types.

/// Simulation Error type.
#[derive(thiserror::Error, Debug)]
pub enum COMError {
    /// Shape Mismatch.
    #[error("Shape Mismatch")]
    ShapeMismatch,

    /// Bad flux coordinate passed to evaluation method.
    #[error("Bad flux coordinate passed to evaluation method")]
    BadFluxCoordinate,

    /// From [`dexter_machine::MachineError`].
    #[error("{0}")]
    MachineError(#[from] dexter_machine::MachineError),

    /// From [`dexter_machine::EvalError`].
    #[error("{0}")]
    EvalError(#[from] dexter_machine::EvalError),
}
