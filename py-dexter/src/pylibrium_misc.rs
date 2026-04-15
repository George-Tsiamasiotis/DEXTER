//! `dexter-equilibrium` newtypes, constructors and method exports.

use crate::{py_debug_impl, py_repr_impl};
use dexter::dexter_equilibrium::LastClosedFluxSurface;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass(name = "_PyLastClosedFluxSurface", frozen)]
pub struct PyLastClosedFluxSurface(pub(crate) LastClosedFluxSurface);

#[pymethods]
impl PyLastClosedFluxSurface {
    #[new]
    pub fn new(kind: &str, value: f64) -> PyResult<Self> {
        Ok(Self(match kind.to_lowercase().as_str() {
            "toroidal" => LastClosedFluxSurface::Toroidal(value),
            "poloidal" => LastClosedFluxSurface::Poloidal(value),
            _ => {
                return Err(PyErr::new::<PyTypeError, _>(
                    "Invalid 'LastClosedFluxSurface'",
                ));
            }
        }))
    }

    #[getter]
    pub fn get_kind(&self) -> &str {
        match self.0 {
            LastClosedFluxSurface::Toroidal(_) => "Toroidal",
            LastClosedFluxSurface::Poloidal(_) => "Poloidal",
        }
    }

    #[getter]
    pub fn get_value(&self) -> f64 {
        self.0.value()
    }
}

py_debug_impl!(PyLastClosedFluxSurface);
py_repr_impl!(PyLastClosedFluxSurface);
