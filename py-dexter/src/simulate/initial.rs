//! Defines wrappers associated with a Particle's and a Queue's initial conditions.

use crate::*;
use dexter::dexter_simulate::*;

use ndarray::Array1;
use numpy::PyReadonlyArray1;
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

#[pyclass(name = "_PyMagneticFluxArray", frozen, immutable_type)]
pub struct PyMagneticFluxArray(Array1<MagneticFlux>);

#[pymethods]
impl PyMagneticFluxArray {
    #[classmethod]
    pub fn toroidal<'py>(_: &Bound<'py, PyType>, array: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let array = array.as_slice()?;
        Ok(Self(toroidal_fluxes(array).into()))
    }

    #[classmethod]
    pub fn poloidal<'py>(_: &Bound<'py, PyType>, array: PyReadonlyArray1<f64>) -> PyResult<Self> {
        let array = array.as_slice()?;
        Ok(Self(poloidal_fluxes(array).into()))
    }
}

// ===============================================================================================

#[pyclass(name = "_PyQueueInitialConditions", frozen, immutable_type)]
pub struct PyQueueInitialConditions(QueueInitialConditions);

#[pymethods]
impl PyQueueInitialConditions {
    #[classmethod]
    pub fn boozer(
        _: &Bound<'_, PyType>,
        t0: PyReadonlyArray1<f64>,
        flux0: &PyMagneticFluxArray,
        theta0: PyReadonlyArray1<f64>,
        zeta0: PyReadonlyArray1<f64>,
        rho0: PyReadonlyArray1<f64>,
        mu0: PyReadonlyArray1<f64>,
    ) -> Result<Self> {
        let flux0_standard_layout = flux0.0.as_standard_layout();
        let flux0_slice = flux0_standard_layout
            .as_slice()
            .expect("probably is in standard layout");
        let initial = QueueInitialConditions::boozer(
            t0.as_slice()?,
            flux0_slice,
            theta0.as_slice()?,
            zeta0.as_slice()?,
            rho0.as_slice()?,
            mu0.as_slice()?,
        )?;
        Ok(Self(initial))
    }

    #[classmethod]
    pub fn mixed(
        _: &Bound<'_, PyType>,
        t0: PyReadonlyArray1<f64>,
        flux0: &PyMagneticFluxArray,
        theta0: PyReadonlyArray1<f64>,
        zeta0: PyReadonlyArray1<f64>,
        pzeta0: PyReadonlyArray1<f64>,
        mu0: PyReadonlyArray1<f64>,
    ) -> Result<Self> {
        let flux0_standard_layout = flux0.0.as_standard_layout();
        let flux0_slice = flux0_standard_layout
            .as_slice()
            .expect("probably is in standard layout");
        let initial = QueueInitialConditions::mixed(
            t0.as_slice()?,
            flux0_slice,
            theta0.as_slice()?,
            zeta0.as_slice()?,
            pzeta0.as_slice()?,
            mu0.as_slice()?,
        )?;
        Ok(Self(initial))
    }
}

// ===============================================================================================

wrapper_debug_export!(PyInitialConditions);
wrapper_debug_export!(PyMagneticFluxArray);
wrapper_debug_export!(PyQueueInitialConditions);

impl_py_repr!(PyInitialConditions, pretty);
impl_py_repr!(PyMagneticFluxArray, pretty);
impl_py_repr!(PyQueueInitialConditions, pretty);
