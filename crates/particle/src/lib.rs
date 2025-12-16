//! # Particle
//!
//! Provides the [`Particle`] object, which represents a Particle inside an [`equilibrium`].
//!
//! The Particle corresponds to a proton with `m=q=1`. All its related quantites and calculated
//! time series are in *Normalized Units*.
//!
//! ## Routines
//!
//! + [`Integration`](Particle::integrate): Straight-forward integration
//! + [`Mapping`](Particle::map): Calculation of the particle's intersections with a θ=const or
//!   ζ=const surface.
//! + [`Single Period Integration`](Particle::single_period_integrate): Straight forward
//!   integration for a single θ-ψp period.
//!
//! ## Integration Configuration
//!
//! Each integration routine's set of parameters can be adjusted with the following structs:
//!
//! + [`IntegrationConfig`]
//! + [`MappingConfig`]
//! + [`SinglePeriodConfig`]

mod config;
mod error;
mod evolution;
mod particle;
mod rkf45;
mod routines;
mod state;

pub(crate) use rkf45::Stepper;

pub use config::*;
pub use error::ParticleError;
pub use evolution::Evolution;
pub use particle::{InitialConditions, IntegrationStatus, OrbitType, Particle};
pub use routines::{Frequencies, MappingParameters, PoincareSection};
pub use state::State;

pub type Result<T> = std::result::Result<T, ParticleError>;
