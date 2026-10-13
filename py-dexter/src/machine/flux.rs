use dexter::dexter_machine::*;
use std::ops::{Add, Div, Mul, Neg, Sub};

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;

// ===============================================================================================

#[pyclass(name = "_PyMagneticFlux", frozen, immutable_type)]
#[derive(PartialEq)]
pub struct PyMagneticFlux(pub(crate) MagneticFlux);

#[pymethods]
impl PyMagneticFlux {
    #[classmethod]
    pub fn toroidal(_: &Bound<'_, PyType>, value: f64) -> Result<Self> {
        Ok(Self(MagneticFlux::Toroidal(value)))
    }

    #[classmethod]
    pub fn poloidal(_: &Bound<'_, PyType>, value: f64) -> Result<Self> {
        Ok(Self(MagneticFlux::Poloidal(value)))
    }

    #[getter]
    pub fn value(&self) -> f64 {
        self.0.value()
    }

    #[getter]
    pub fn kind(&self) -> String {
        format!("{:?}", self.0.kind())
    }

    // Unfortunately, `to_other` cannot be implemented as it leads to circular imports in python
}

#[pymethods]
impl PyMagneticFlux {
    fn __add__(&self, scalar: f64) -> PyMagneticFlux {
        PyMagneticFlux(self.0.add(scalar))
    }

    fn __sub__(&self, scalar: f64) -> PyMagneticFlux {
        PyMagneticFlux(self.0.sub(scalar))
    }

    fn __neg__(&self) -> PyMagneticFlux {
        PyMagneticFlux(self.0.neg())
    }

    fn __mul__(&self, scalar: f64) -> PyMagneticFlux {
        PyMagneticFlux(self.0.mul(scalar))
    }

    fn __truediv__(&self, scalar: f64) -> PyMagneticFlux {
        PyMagneticFlux(self.0.div(scalar))
    }
}

impl From<MagneticFlux> for PyMagneticFlux {
    fn from(value: MagneticFlux) -> Self {
        Self(value)
    }
}

// ===============================================================================================

#[pyclass(name = "_PyMagneticFluxArray", immutable_type, frozen)]
pub struct PyMagneticFluxArray(pub MagneticFluxArray);

#[pymethods]
impl PyMagneticFluxArray {
    #[classmethod]
    pub fn toroidal<'py>(_: &Bound<'py, PyType>, array: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let array = array.as_array().as_standard_layout().to_owned();
        Ok(Self(MagneticFluxArray::Toroidal(array)))
    }

    #[classmethod]
    pub fn poloidal<'py>(_: &Bound<'py, PyType>, array: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let array = array.as_array().as_standard_layout().to_owned();
        Ok(Self(MagneticFluxArray::Poloidal(array)))
    }

    #[getter]
    pub fn values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.values().to_owned().into_pyarray(py)
    }

    #[getter]
    pub fn kind(&self) -> String {
        format!("{:?}", self.0.kind())
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }

    pub fn __getitem__(&self, index: usize) -> Option<PyMagneticFlux> {
        self.0.get(index).map(|flux| PyMagneticFlux(flux))
    }

    // Unfortunately, `to_other` cannot be implemented as it leads to circular imports in python
}

// ===============================================================================================

wrapper_debug_export!(PyMagneticFlux);
wrapper_debug_export!(PyMagneticFluxArray);
impl_py_repr!(PyMagneticFlux, simple);
impl_py_repr!(PyMagneticFluxArray, pretty);
