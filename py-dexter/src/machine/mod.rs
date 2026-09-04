//! Defines `dexter-machine` helper objects and exports machine objects and helper types

mod bfield;
mod current;
mod geometry;
mod mode;
mod perturbation;
mod qfactor;

pub use bfield::*;
pub use current::*;
pub use geometry::*;
pub use mode::*;
pub use perturbation::*;
pub use qfactor::*;

use crate::{impl_py_repr, wrapper_debug_export};
use dexter::dexter_machine::*;
use pyo3::{prelude::*, types::PyType};

// ===============================================================================================

#[pyclass(name = "_PyMagneticFlux", frozen, immutable_type)]
#[derive(PartialEq)]
pub struct PyMagneticFlux(pub(crate) MagneticFlux);

#[pymethods]
impl PyMagneticFlux {
    #[classmethod]
    pub fn toroidal(_: &Bound<'_, PyType>, value: f64) -> PyResult<Self> {
        Ok(Self(MagneticFlux::Toroidal(value)))
    }

    #[classmethod]
    pub fn poloidal(_: &Bound<'_, PyType>, value: f64) -> PyResult<Self> {
        Ok(Self(MagneticFlux::Poloidal(value)))
    }

    #[getter]
    pub fn value(&self) -> f64 {
        self.0.value()
    }

    #[getter]
    pub fn kind(&self) -> String {
        self.0.kind().into()
    }
}

impl From<MagneticFlux> for PyMagneticFlux {
    fn from(value: MagneticFlux) -> Self {
        Self(value)
    }
}

wrapper_debug_export!(PyMagneticFlux);
impl_py_repr!(PyMagneticFlux, simple);
