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
    pub fn eval_ms(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_ms(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dms_dflux(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dms_dflux(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dms_dtheta(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dms_dtheta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dms_dzeta(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dms_dzeta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dms_dt(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dms_dt(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_phis(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_phis(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dphis_dflux(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dphis_dflux(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dphis_dtheta(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dphis_dtheta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dphis_dzeta(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dphis_dzeta(flux, theta, zeta, t, &mut self.0.generate_caches())?)
    }

    pub fn eval_dphis_dt(&self, theta: f64, zeta: f64, t: f64, psit:f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.0.eval_dphis_dt(flux, theta, zeta, t, &mut self.0.generate_caches())?)
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
