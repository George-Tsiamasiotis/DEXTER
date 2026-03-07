//! Simulations of charged particles inside an [`Equilibrium`](dexter_equilibrium)
//!
//! ### [`Particle`]: A charged particle inside an [`Equilibrium`](dexter_equilibrium)
//!
//! The Particle corresponds to a proton with `m=q=1`. All its related quantities and calculated
//! time series are in *Normalized Units*.
//!
//! Routines:
//!
//! + [`Particle::integrate`]: Integrates the particle in a specific time interval.
//!
//! #### Integration Configuration
//!
//! Each integration routine's set of parameters can be adjusted with the following structs:
//!
//! + [`SolverParams`]: Parameters passed to the solver.
//!
//! The integration is done with the RKF4(5) method. The step size can be adaptive by minimizing
//! the energy difference from step to step or minimizing the local truncation error (classic
//! RKF4(5)), or simple set to be constant.

#![allow(confusable_idents)]

mod error;
mod particle;
mod state;
mod system;

// ============== Public API

pub use error::SimulationError;
/// Result type used throughout the crate.
pub type Result<T> = std::result::Result<T, SimulationError>;

pub(crate) use system::FluxCoordinate;
pub use system::{SolverParams, SteppingMethod};

pub use particle::{InitialConditions, InitialFlux, IntegrationStatus, Particle};
