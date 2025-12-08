#![doc = include_str!("../README.md")]

mod bfields;
mod cache;
mod currents;
mod error;
mod eval;
mod geometry;
mod harmonics;
mod perturbation;
mod qfactors;

pub mod extract;

pub use eval::{Bfield, Current, Geometry, Harmonic, Perturbation, Qfactor};

pub use bfields::{NcBfield, NcBfieldBuilder};
pub use currents::{NcCurrent, NcCurrentBuilder};
pub use geometry::{NcGeometry, NcGeometryBuilder};
pub use harmonics::{NcHarmonic, NcHarmonicBuilder};
pub use qfactors::{NcQfactor, NcQfactorBuilder};

pub use harmonics::PhaseMethod;
pub use perturbation::NcPerturbation;

pub use cache::{HarmonicCache, NcHarmonicCache};

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

/// Creates a getter method for extracting the flat Vec data as an Array2.
#[doc(hidden)]
#[macro_export]
macro_rules! fortran_vec_to_carray2d_impl {
    ($meth_name:ident, $field:ident, $var_name:ident) => {
        /// Returns the `R(ψp, θ)` data as a 2D array.
        #[doc = "Returns the `"]
        #[doc = stringify!($var_name)]
        #[doc = "`(ψp, θ) data as an [`Array2<f64>`]"]
        pub fn $meth_name(&self) -> Array2<f64> {
            // Array is in Fortran order.
            let shape = (self.psip_data.len(), self.theta_data.len());
            Array2::from_shape_vec(shape, self.$field.to_vec())
                .expect("shape is correct by definition")
                .reversed_axes()
        }
    };
}
