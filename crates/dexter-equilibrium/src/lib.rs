#![doc = include_str!("../README.md")]
//!
//! # Equilibrium objects
//!
//! + Representations of an equilibrium's general geometry. Provides interpolation methods between `ψ`, `ψp`, `r`, `R`, `Z`, `J`.
//!     - [`LarGeometry`]: Analytical Large Aspect Ratio Geometry of a circular device.
//!     - [`NcGeometry`]: Geometry of the netCDF equilibrium
//!
//! + Representations of the q-factor profile:
//!     - [`UnityQfactor`]: q-factor profile of q = 1 and ψ=ψp.
//!     - [`ParabolicQfactor`]: q-factor of parabolic q(ψ) profile.
//!     - [`NcQfactor`]: q-factor reconstructed from a netCDF file.
//!
//! + Representations of the plasma currents:
//!     - [`LarCurrent`]: Large Aspect Ration plasma current with g=1 and I=0.
//!     - [`NcCurrent`]: Plasma current reconstructed from a netCDF file.
//!
//! + Representations of the magnetic field:
//!     - [`LarBfield`]: Large Aspect Ratio magnetic field with B(ψ, θ) = 1 - sqrt(2ψ)cos(θ).
//!     - [`NcBfield`]: Magnetic reconstructed from a netCDF file.
//!
//! + Representations of single perturbation harmonics:
//!     - [`CosHarmonic`]: Simple analytical harmonic of the form `α*cos(mθ-nζ+φ)`.
//!     - [`NcHarmonic`]: Single perturbation harmonic from a netCDF file of the form
//!     `α(ψ/ψp) * cos(mθ-nζ+φ(ψ/ψp))`
//!
//! + Representations of Perturbations.
//!     - [`Perturbation`]: A sum of an arbitrary number of [`Harmonics`](Harmonic).
//!
//! ## Evaluations:
//!
//! + [`Geometry`]: Conversions to laboratory quantities.
//! + [`FluxCommute`]: Conversion between the two flux coordinates `ψ` and `ψp`
//! + [`Qfactor`]: Evaluation of q-factor related quantities.
//! + [`Current`]: Evaluation of plasma current related quantities.
//! + [`Bfield`]: Evaluation of magnetic field related quantities.
//! + [`Harmonic`]: Single perturbation harmonic related quantities computation.
//!
//! ## Caching
//!
//! The trait [`Harmonic`] requires a [`HarmonicCache`] object to be passed as a parameter. Such an
//! object caches values such as angles' modulos and their sines/cosines, or amplitudes/phases
//! calculated with interpolation. It may also store the necessary
//! [`Accelerators`](rsl_interpolation::Accelerator). Since many evaluation methods are called with
//! the same arguments sequentially (like in the case of the perturbed equations of motion), it
//! makes sense to cache values that appear many times in these methods and can be expensive in
//! tight loops.
//!
//! + [`CosHarmonicCache`]: Cache for [`CosHarmonic`]
//! + [`NcHarmonicCache`]: Cache for [`NcHarmonic`]
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

#![allow(mixed_script_confusables)]

mod error;
mod eval;
mod flux;
mod objects;

// ============== Public API

pub mod extract;

pub use error::EqError;
pub use error::NcError;
/// Result type used throughout the crate.
pub type Result<T> = std::result::Result<T, EqError>;

pub use objects::EquilibriumType;

pub use flux::NcFlux;
pub use flux::NcFluxState;

pub use eval::HarmonicCache;
pub use eval::{Bfield, Current, FluxCommute, Geometry, Harmonic, Qfactor};

pub use objects::geometries::LarGeometry;
pub use objects::geometries::NcGeometry;
pub use objects::geometries::NcGeometryBuilder;

pub use objects::qfactors::NcQfactor;
pub use objects::qfactors::NcQfactorBuilder;
pub use objects::qfactors::UnityQfactor;
pub use objects::qfactors::{FluxWall, ParabolicQfactor};

pub use objects::currents::LarCurrent;
pub use objects::currents::NcCurrent;
pub use objects::currents::NcCurrentBuilder;

pub use objects::bfield::LarBfield;
pub use objects::bfield::NcBfield;
pub use objects::bfield::NcBfieldBuilder;

pub use objects::harmonics::{CosHarmonic, CosHarmonicCache};
pub use objects::nc_harmonic::{NcHarmonic, NcHarmonicBuilder, NcHarmonicCache, PhaseMethod};

pub use objects::perturbation::Perturbation;
