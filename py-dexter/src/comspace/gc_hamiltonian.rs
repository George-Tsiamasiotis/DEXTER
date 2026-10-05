//! Defines the `GCHamiltonian` wrapper.

use crate::*;
use dexter::dexter_comspace::*;
use dexter::dexter_machine::*;

use numpy::PyReadonlyArray1;
use numpy::PyReadonlyArray2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyType;

#[pyclass(name = "_PyGCHamiltonian", immutable_type, frozen)]
pub struct PyGCHamiltonian;

// ===============================================================================================

#[pymethods]
impl PyGCHamiltonian {
    #[new]
    pub fn new() -> Self {
        Self
    }

    #[classmethod]
    pub fn eval(
        _: &Bound<'_, PyType>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta: f64,
        flux: &PyMagneticFlux,
        pzeta: f64,
        mu: f64,
    ) -> Result<f64> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        Ok(GCHamiltonian::eval(machine, theta, flux.0, pzeta, mu)?)
    }

    #[classmethod]
    pub fn eval_psi_array<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_array: PyReadonlyArray1<f64>,
        psi_array: PyReadonlyArray1<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray1<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let grid = GCHamiltonian::eval_psi_array(
            machine,
            &theta_array.as_array(),
            &psi_array.as_array(),
            pzeta,
            mu,
        )?;
        Ok(grid.into_pyarray(py))
    }

    #[classmethod]
    pub fn eval_psip_array<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_array: PyReadonlyArray1<f64>,
        psip_array: PyReadonlyArray1<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray1<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let grid = GCHamiltonian::eval_psip_array(
            machine,
            &theta_array.as_array(),
            &psip_array.as_array(),
            pzeta,
            mu,
        )?;
        Ok(grid.into_pyarray(py))
    }

    #[classmethod]
    pub fn eval_psi_grid<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_grid: PyReadonlyArray2<f64>,
        psi_grid: PyReadonlyArray2<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let grid = GCHamiltonian::eval_psi_grid(
            machine,
            &theta_grid.as_array(),
            &psi_grid.as_array(),
            pzeta,
            mu,
        )?;
        Ok(grid.into_pyarray(py))
    }

    #[classmethod]
    pub fn eval_psip_grid<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_grid: PyReadonlyArray2<f64>,
        psip_grid: PyReadonlyArray2<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let grid = GCHamiltonian::eval_psip_grid(
            machine,
            &theta_grid.as_array(),
            &psip_grid.as_array(),
            pzeta,
            mu,
        )?;
        Ok(grid.into_pyarray(py))
    }
}

// ===============================================================================================

#[pyfunction]
#[pyo3(name = "_py_create_poloidal_grid")]
pub fn py_create_poloidal_grid<'py>(
    py: Python<'py>,
    theta_array: PyReadonlyArray1<f64>,
    flux_array: PyReadonlyArray1<f64>,
) -> (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>) {
    let theta_array = theta_array.as_array();
    let flux_array = flux_array.as_array();
    let (theta_grid, flux_grid) = create_poloidal_grid(&theta_array, &flux_array);

    (
        theta_grid.to_owned().into_pyarray(py),
        flux_grid.to_owned().into_pyarray(py),
    )
}

// ===============================================================================================

impl std::fmt::Debug for PyGCHamiltonian {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("GCHamiltonian")
    }
}

#[pymethods]
impl PyGCHamiltonian {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self)
    }
}
