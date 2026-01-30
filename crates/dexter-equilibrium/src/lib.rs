#![doc = include_str!("../README.md")]
//!
//! # Equilibrium objects
//!
//! + [`Current`]: Representation of the plasma currents
//!     - [`NcCurrent`]: Plasma current reconstructed from a netCDF file.
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
//!
//! ## Caching
//!
//! TODO:

mod error;
mod eval;
mod flux;
mod objects;

pub(crate) use flux::*;

pub mod extract;

pub use eval::Current;
pub use objects::*;

pub use error::{EqError, NcError};
pub type Result<T> = std::result::Result<T, EqError>;
