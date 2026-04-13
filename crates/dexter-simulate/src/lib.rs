//! Simulations of charged particles inside an [`Equilibrium`](dexter_equilibrium).
//!
//! ### [`Particle`]: A charged particle inside an [`Equilibrium`](dexter_equilibrium)
//!
//! The Particle corresponds to a proton with `m=q=1`. All its related quantities and calculated
//! time series are in *Normalized Units*.
//!
//! Routines:
//!
//! + [`Particle::integrate`]: Integrates the particle in a specific time interval.
//! + [`Particle::intersect`]: Integrates the particle, calculating its intersections with a
//!   constant `θ` or `ζ` surface.
//!
//! ### [`Queue`]: A container for batch initializing and running Particle routines.
//!
//! A Queue is conveniently initialized from 1D [`arrays`](ndarray::ArrayBase) of initial conditions with
//! the help of the [`QueueInitialConditions`] helper struct.
//!
//! Routines:
//!
//! + [`Queue::integrate`]: Integrates the contained particles in a specific time interval.
//! + [`Queue::intersect`]: Integrates the contained particles, calculating their intersections with a
//!   constant `θ` or `ζ` surface. Otherwise known as a `Poincare map`.
//!
//! ### Parallelism
//!
//! Using the [`rayon`] crate, particle routines can run in parallel. The number of threads to be
//! used can be adjusted with the [`set_num_threads`] function.
//!
//! ### Integration Configuration
//!
//! Each integration routine's set of parameters can be adjusted with the following structs:
//!
//! + [`SolverParams`]: Parameters passed to the solver.
//!
//! The integration is done with the RKF4(5) method. The step size can be adaptive by minimizing
//! the energy difference from step to step or minimizing the local truncation error (classic
//! RKF4(5)), or simple set to be constant.
//!
//! ### Constants of motion
//!
//! Calculations involving the Constants of Motion.
//!
//! + [`COMs`]: A container struct for the constants of motion in an unperturbed equilibrium.
//!   + [`COMs::energy_of_psi_grid`]: Calculation of the Energy in a 2D `θ-ψ` grid
//!   + [`COMs::energy_of_psip_grid`]: Calculation of the Energy in a 2D `θ-ψp` grid
//!

mod coms;
mod error;
mod particle;
mod queue;
mod solve;
mod state;

// ============== Public API

pub use dexter_common::{get_max_threads, set_num_threads};

pub use error::SimulationError;
/// Result type used throughout the crate.
pub type Result<T> = std::result::Result<T, SimulationError>;

pub use solve::{FluxCoordinate, SolverParams, SteppingMethod};

pub use particle::{
    CoordinateSet, Frequencies, InitialConditions, InitialFlux, IntegrationStatus, IntersectParams,
    Intersection, OrbitType, Particle, ParticleCacheStats,
};

pub use queue::{Queue, QueueInitialConditions, Routine};
pub use queue::{poloidal_fluxes, toroidal_fluxes};

pub use coms::COMs;

// ============== Configuration constants

/// Crate configuration constants.
pub mod constants {
    use super::SteppingMethod;

    /// The default optimal step calculation method.
    pub const DEFAULT_STEPPING_METHOD: SteppingMethod = SteppingMethod::EnergyAdaptiveStep;

    /// The default maximum amount of steps a particle can make before terminating its integration.
    pub const DEFAULT_MAX_STEPS: usize = 1_000_000;

    /// The default initial time step for the RKF45 adaptive step method. The value is empirical.
    pub const DEFAULT_FIRST_STEP: f64 = 1e-1;

    /// The default safety factor of the solver. Should be less than 1.0.
    pub const DEFAULT_SAFETY_FACTOR: f64 = 0.9;

    /// The default relative tolerance of the energy difference in every step.
    pub const DEFAULT_ENERGY_REL_TOL: f64 = 1e-12;

    /// The default absolute tolerance of the energy difference in every step.
    pub const DEFAULT_ENERGY_ABS_TOL: f64 = 1e-14;

    /// The default relative tolerance of the local truncation error in every step.
    pub const DEFAULT_ERROR_REL_TOL: f64 = 1e-12;

    /// The default absolute tolerance of the local truncation error in every step.
    pub const DEFAULT_ERROR_ABS_TOL: f64 = 1e-14;

    /// The default relative tolerance of the `ψ-ψ0` check in the
    /// [`Particle::close`][crate::Particle::close] routine.
    pub const PSI_RELATIVE_TOLERANCE: f64 = 1e-6;

    /// See [`OrbitType::TrappedStagnated`][crate::OrbitType::TrappedStagnated].
    pub const TRAPPED_STAGNATED_CLASSIFICATION_CUTOFF: f64 = 0.999;
}
