/// Custom error types
#[derive(thiserror::Error, Debug)]
pub enum ParticleError {
    /// From [`equilibrium::EqError`].
    #[error("{0}")]
    // FIXME: At the moment, this can only mean DomainError, but if other cases appear in the
    // future, we must distinguish between the different error variants.
    EqError(#[from] equilibrium::EqError),

    /// Particle timed out.
    #[error("Particle timed out after {0:?}")]
    TimedOut(std::time::Duration),

    /// Intersection accuracy check failed.
    #[error("Intersection accuracy check failed.")]
    IntersectionError,

    /// NaN values encountered in State evaluation.
    #[error("NaN values encountered in State evaluation.")]
    EvaluationNaN,
}
