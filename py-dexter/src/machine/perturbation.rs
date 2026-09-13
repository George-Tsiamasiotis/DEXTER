//! Defines the PyPerturbation container.

use std::sync::Arc;

use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(
    name = "_PyPerturbation",
    frozen,
    immutable_type,
    sequence,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPerturbation(Arc<Perturbation>);

#[pymethods]
impl PyPerturbation {
    #[new]
    pub fn new<'py>(modes: Bound<'py, PyList>) -> Result<Self> {
        let pymodes: Vec<PyMode> = modes.into_iter().map(|m| m.extract().unwrap()).collect();
        let modes: DynModes = pymodes.iter().map(|m| m.boxed_mode()).collect();
        let perturbation = Perturbation::new(modes);
        Ok(PyPerturbation(Arc::new(perturbation)))
    }

    pub fn __len__(&self) -> usize {
        self.0.count()
    }
}

impl PyPerturbation {
    pub fn inner(&self) -> &Perturbation {
        &self.0
    }
}

// ===============================================================================================

#[pymethods] // Evaluations
#[rustfmt::skip]
impl PyPerturbation {
    pub fn eval_p(&self, theta: f64, zeta: f64, t: f64, psi:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.0.eval_p(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_deriv_flux(&self, theta: f64, zeta: f64, t: f64, psi:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.0.eval_deriv_flux(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_deriv_theta(&self, theta: f64, zeta: f64, t: f64, psi:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.0.eval_deriv_theta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_deriv_zeta(&self, theta: f64, zeta: f64, t: f64, psi:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.0.eval_deriv_zeta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_deriv_t(&self, theta: f64, zeta: f64, t: f64, psi:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.0.eval_deriv_t(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }
}

// ===============================================================================================

wrapper_debug_export!(PyPerturbation);

#[pymethods]
impl PyPerturbation {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.0)
    }
}
