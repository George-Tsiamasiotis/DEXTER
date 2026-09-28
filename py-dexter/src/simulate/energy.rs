use dexter::dexter_machine::MachineBuilder;
use dexter::dexter_simulate;
use numpy::{IntoPyArray, PyReadonlyArray1, PyReadonlyArray2};

use crate::*;
use pyo3::prelude::*;

use numpy::PyArray2;

#[pyfunction]
#[pyo3(name = "_py_create_poloidal_grid")]
pub fn create_poloidal_grid<'py>(
    py: Python<'py>,
    theta_array: PyReadonlyArray1<f64>,
    flux_array: PyReadonlyArray1<f64>,
) -> (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>) {
    let theta_array = theta_array.as_array();
    let flux_array = flux_array.as_array();
    let (theta_grid, flux_grid) = dexter_simulate::create_poloidal_grid(&theta_array, &flux_array);

    (
        theta_grid.to_owned().into_pyarray(py),
        flux_grid.to_owned().into_pyarray(py),
    )
}

#[pyfunction]
#[pyo3(name = "_py_energy_of_psi_grid")]
pub fn energy_of_psi_grid<'py>(
    py: Python<'py>,
    qfactor: &PyQfactor,
    current: &PyCurrent,
    bfield: &PyBfield,
    perturbation: &PyPerturbation,
    pzeta: f64,
    mu: f64,
    theta_array: PyReadonlyArray2<f64>,
    psi_array: PyReadonlyArray2<f64>,
) -> Bound<'py, PyArray2<f64>> {
    let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
        .with_perturbation(perturbation.inner())
        .build();
    let energy_grid = dexter_simulate::energy_of_psi_grid(
        machine,
        pzeta,
        mu,
        &theta_array.as_array(),
        &psi_array.as_array(),
    );
    energy_grid.into_pyarray(py)
}

#[pyfunction]
#[pyo3(name = "_py_energy_of_psip_grid")]
pub fn energy_of_psip_grid<'py>(
    py: Python<'py>,
    qfactor: &PyQfactor,
    current: &PyCurrent,
    bfield: &PyBfield,
    perturbation: &PyPerturbation,
    pzeta: f64,
    mu: f64,
    theta_array: PyReadonlyArray2<f64>,
    psip_array: PyReadonlyArray2<f64>,
) -> Bound<'py, PyArray2<f64>> {
    let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
        .with_perturbation(perturbation.inner())
        .build();
    let energy_grid = dexter_simulate::energy_of_psip_grid(
        machine,
        pzeta,
        mu,
        &theta_array.as_array(),
        &psip_array.as_array(),
    );
    energy_grid.into_pyarray(py)
}
