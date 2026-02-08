#![doc = include_str!("../README.md")]
//!
//! # Equilibrium objects
//!
//! + Representations of an equilibrium's general geometry. Provides interpolation methods between `ψ`, `ψp`, `r`, `R`, `Z`, `J`.
//!     - [`NcGeometry`]: Geometry of the netCDF equilibrium
//!
//! + Representations of the q-factor profile:
//!     - [`NcQfactor`]: q-factor reconstructed from a netCDF file.
//!     - [`UnityQfactor`]: q-factor profile of q = 1 and ψ=ψp.
//!     - [`ParabolicQfactor`]: q-factor of parabolic q(ψ) profile.
//!
//! + Representations of the plasma currents;
//!     - [`NcCurrent`]: Plasma current reconstructed from a netCDF file.
//!     - [`LarCurrent`]: Large Aspect Ration plasma current with g=1 and I=0.
//!
//! ## Evaluations:
//!
//! + [`Geometry`]: Conversions to laboratory quantities.
//! + [`FluxCommute`]: Conversion between the two flux coordinates `ψ` and `ψp`
//! + [`Qfactor`]: Evaluation of q-factor related quantities.
//! + [`Current`]: Evaluation of plasma current related quantities.
//!
//! TODO:
//!
//! ## Data extraction
//!
//! The [`extract`] module provides methods for extracting data arrays and
//! [`Variables`](netcdf::Variable) from the netCDF file.
//!
//! + [`extract::open`]: Open a netCDF file.
//! + [`extract::extract_scalar`]: Scalar value extraction.
//! + [`extract::extract_1d_array`]: 1D array extraction.
//! + [`extract::extract_2d_array`]: 2D array extraction.
//! + [`extract::extract_3d_array`]: 3D array extraction.
//! + [`extract::extract_harmonic_arrays`]: Extraction of the α and φ arrays of the {m,n} mode.
//! Modes are indexed by their mode numbers, rather than the logical index they appear on the data
//! arrays.
//! + [`extract::extract_variable`]: Extraction of a variable as a [`Variable`](netcdf::Variable).
//! + [`extract::extract_attribute`]: Extraction of a file's attribute as a String.
//! + [`extract::extract_version`]: Extraction of a files convention [`Semantic Version`](https://semver.org/)
//!
//! ## Caching
//!
//! TODO:

#![allow(mixed_script_confusables)]

mod error;
mod eval;
mod flux;
mod objects;

// ============== Public API

pub mod extract;

pub use error::EqError;
pub use error::NcError;
pub type Result<T> = std::result::Result<T, EqError>;

pub use objects::EquilibriumType;

pub use flux::NcFlux;
pub use flux::NcFluxState;

pub use eval::{Current, FluxCommute, Geometry, Qfactor};

pub use objects::geometries::NcGeometry;
pub use objects::geometries::NcGeometryBuilder;

pub use objects::qfactors::NcQfactor;
pub use objects::qfactors::NcQfactorBuilder;
pub use objects::qfactors::UnityQfactor;
pub use objects::qfactors::{FluxWall, ParabolicQfactor};

pub use objects::currents::LarCurrent;
pub use objects::currents::NcCurrent;
pub use objects::currents::NcCurrentBuilder;
