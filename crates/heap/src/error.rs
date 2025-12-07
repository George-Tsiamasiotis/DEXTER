/// Custom error types for a heap of Particles
#[derive(thiserror::Error, Debug)]
pub enum HeapError {
    /// From [`particle::ParticleError`].
    #[error("{0}")]
    ParticleError(#[from] particle::ParticleError),

    /// ndarray's ShapeError.
    #[error("Shape Error: {0}")]
    ShapeError(#[from] ndarray::ShapeError),

    /// Heap initial conditions arrays must be of the same size.
    #[error("Initial conditions arrays must be of the same size")]
    InitMismatch,
}
