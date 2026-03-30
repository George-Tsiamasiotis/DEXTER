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
//! Using the [`rayon`] crate, particle routines can run in parallel.
//!
//! Routines:
//!
//! + [`Queue::integrate`]: Integrates the contained particles in a specific time interval.
//! + [`Queue::intersect`]: Integrates the contained particles, calculating their intersections with a
//!   constant `θ` or `ζ` surface. Otherwise known as a `Poincare map`.
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

mod error;
mod particle;
mod queue;
mod solve;
mod state;

// ============== Public API

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

pub use state::COMs;
