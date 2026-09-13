//! Defines wrappers associated with a Particle's and a Queue's initial conditions.

use crate::*;
use dexter::dexter_simulate::*;

use pyo3::{prelude::*, types::PyType};

// ===============================================================================================

#[pyclass(name = "_PyInitialConditions", frozen, immutable_type)]
pub struct PyInitialConditions(pub(crate) InitialConditions);

#[pymethods]
impl PyInitialConditions {
    #[classmethod]
    pub fn boozer(
        _: &Bound<'_, PyType>,
        t0: f64,
        flux0: &PyMagneticFlux,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> Self {
        Self(InitialConditions::boozer(
            t0, flux0.0, theta0, zeta0, rho0, mu0,
        ))
    }

    #[classmethod]
    pub fn mixed(
        _: &Bound<'_, PyType>,
        t0: f64,
        flux0: &PyMagneticFlux,
        theta0: f64,
        zeta0: f64,
        pzeta0: f64,
        mu0: f64,
    ) -> Self {
        Self(InitialConditions::mixed(
            t0, flux0.0, theta0, zeta0, pzeta0, mu0,
        ))
    }

    #[getter]
    pub fn t0(&self) -> f64 {
        self.0.t0()
    }

    #[getter]
    pub fn flux0(&self) -> PyMagneticFlux {
        self.0.flux0().into()
    }

    #[getter]
    pub fn theta0(&self) -> f64 {
        self.0.theta0()
    }

    #[getter]
    pub fn zeta0(&self) -> f64 {
        self.0.zeta0()
    }

    #[getter]
    pub fn rho0(&self) -> Option<f64> {
        self.0.rho0()
    }

    #[getter]
    pub fn pzeta0(&self) -> Option<f64> {
        self.0.pzeta0()
    }

    #[getter]
    pub fn mu0(&self) -> f64 {
        self.0.mu0()
    }

    #[getter]
    pub fn coordinate_set(&self) -> String {
        format!("{:?}", self.0.coordinate_set())
    }
}

// ===============================================================================================

wrapper_debug_export!(PyInitialConditions);

#[pymethods]
impl PyInitialConditions {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self)
    }
}
