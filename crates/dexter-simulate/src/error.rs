/// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum SimulationError {
    /// From [`dexter_equilibrium::EqError`].
    #[error("{0}")]
    EqError(#[from] dexter_equilibrium::EqError),

    /// NaN encountered inside the Stepper.
    #[error("NaN encountered inside the Stepper")]
    StepperNan,
}
