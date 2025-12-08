#![doc = include_str!("../README.md")]

mod cache;
mod error;
mod eval;
pub mod extract;
mod objects;

pub use objects::*;

pub use cache::{HarmonicCache, NcHarmonicCache};
pub use eval::{Bfield, Current, Geometry, Harmonic, Perturbation, Qfactor};

pub use error::{EqError, NcError};

pub type Result<T> = std::result::Result<T, EqError>;

/// Magnetic flux, in Normalized Units.
#[doc(alias = "f64")]
pub type Flux = f64;

/// Angle in radians.
#[doc(alias = "f64")]
pub type Radians = f64;

/// Distance, in Normalized Units (normalized to the major radius R).
#[doc(alias = "f64")]
pub type Length = f64;
