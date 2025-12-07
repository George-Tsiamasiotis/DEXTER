#![doc = include_str!("../README.md")]
//!
//! # Equilibrium objects
//!
//! + [`Qfactor`]: Representation of the q-factor profile
//!     - [`qfactors::NcQfactor`]: q-factor reconstructed from a netCDF file.
//!     - [`qfactors::Unity`]: Analytical q-factor profile with q=1.
//!
//! + [`Current`]: Representation of the plasma currents
//!     - [`currents::NcCurrent`]: Plasma current reconstructed from a netCDF file.
//!     - [`currents::LarCurrent`]: Large Aspect Ratio approximation, with g=1 and I=0.
//!
//! + [`Bfield`]: Representation of the magnetic field
//!     - [`bfields::NcBfield`]: Magnetic field reconstructed from a netCDF file.
//!
//! + [`Harmonic`]: Representation of a single magnetic ripple harmonic
//!     - [`harmonics::NcHarmonic`]: Harmonic reconstructed from a netCDF file.
//!         - [`harmonics::PhaseMethod`]: The method for calculating the reconstructed Harmonic's
//!         phase `φ(ψp)`.
//!
//! + [`Perturbation`]: Representation of a sum of different harmonics.
//!     - [`perturbations::NcPerturbation`]: Sum of arbitrarily many [`harmonics::NcHarmonic`]s.
//!
//! ## Geometry object
//!
//! + [`Geometry`]: An object containing all information about the geometry of the equilibrium, and
//! provides interpolation methods between `ψp`, `r`, `R` and `Z`
//!     - [`geometries::NcGeometry`]: Geometry of the netCDF equilibrium
//!
//! ## Data extraction
//!
//! Methods for extracting data arrays and [`Variables`](netcdf::Variable) from the netCDF file.
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
//!
//! ## Caching
//!
//! + [`cache::HarmonicCache`]: Cache holding a single Harmonic's intermediate values.
//!
//! Since each harmonic appears 3-4 times in the equation of motion in the form of its different
//! derivatives, we can cache its *angular* part, as well as its *sin* and *cos* values, to avoid
//! recalculating them at every call. This is an important performance boost, since the sin/cos
//! functions are quite expensive when evaluated so many times in such a tight loop. We can also
//! cache all its interpolated values.

mod error;
mod eval;
mod objects;

pub mod extract;
pub use objects::*;
pub mod cache;

pub use eval::{Bfield, Current, Geometry, Harmonic, Perturbation, Qfactor};

pub use error::{EqError, NcError};

pub type Result<T> = std::result::Result<T, EqError>;
