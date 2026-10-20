//! Exports [`EnergyPzetaPlane`] and its helper types.

use dexter::dexter_comspace::*;
use dexter::dexter_machine::*;
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use parabola::{Line, LineIntercepts, Parabola};
use pyo3::types::PyList;

use crate::*;
use pyo3::prelude::*;

// ===============================================================================================

/// Only used internally
#[pyclass(name = "_PyParabola", immutable_type, frozen)]
pub struct PyParabola(pub Parabola);

#[pymethods]
impl PyParabola {
    pub fn eval_array<'py>(
        &self,
        py: Python<'py>,
        array: PyReadonlyArray1<f64>,
    ) -> Bound<'py, PyArray1<f64>> {
        array.as_array().map(|x| self.0.eval(*x)).into_pyarray(py)
    }

    pub fn horizontal_intercepts(&self, y: f64) -> (f64, f64) {
        let line = Line {
            slope: 0.0,
            intercept: y,
        };
        match self.0.line_intercepts(&line) {
            LineIntercepts::TwoIntercepts(first, second) => (first.x, second.x),
            _ => panic!("Line does not intercept the parabola at 2 points"),
        }
    }
}

// ===============================================================================================

#[pyclass(name = "_PyEnergyPzetaPlane")]
pub struct PyEnergyPzetaPlane(pub EnergyPzetaPlane);

#[pymethods]
impl PyEnergyPzetaPlane {
    #[new]
    pub fn new(qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield, mu: f64) -> Self {
        // EnergyPzetaPlane does not use a Geometry or a Perturbation
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        Self(EnergyPzetaPlane::from_mu(machine, mu))
    }

    pub fn calculate_critical_points(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        pzeta_array: PyReadonlyArray1<f64>,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        self.0
            .calculate_critical_points(machine, &pzeta_array.as_array());
    }

    pub fn calculate_boundary_extrema(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        pzeta_array: PyReadonlyArray1<f64>,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        self.0
            .calculate_boundary_extrema(machine, &pzeta_array.as_array());
    }

    #[getter]
    pub fn mu(&self) -> f64 {
        self.0.mu()
    }

    #[getter]
    pub fn axis_parabola(&self) -> PyParabola {
        PyParabola(self.0.axis_parabola().clone())
    }

    #[getter]
    pub fn left_wall_parabola(&self) -> PyParabola {
        PyParabola(self.0.left_wall_parabola().clone())
    }

    #[getter]
    pub fn right_wall_parabola(&self) -> PyParabola {
        PyParabola(self.0.right_wall_parabola().clone())
    }

    #[getter]
    pub fn tp_pzeta_values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.tp_pzeta_values().to_owned().into_pyarray(py)
    }

    #[getter]
    pub fn tp_upper_values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.tp_upper_values().to_owned().into_pyarray(py)
    }

    #[getter]
    pub fn tp_lower_values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.tp_lower_values().to_owned().into_pyarray(py)
    }

    #[getter]
    pub fn critical_points<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let iter = self
            .0
            .critical_points()
            .iter()
            .map(|critical_point| PyCriticalPoint(critical_point.clone()));
        PyList::new(py, iter)
    }

    #[getter]
    pub fn boundary_extrema<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let iter = self
            .0
            .boundary_extrema()
            .iter()
            .map(|boundary_extremum| PyBoundaryExtremum(boundary_extremum.clone()));
        PyList::new(py, iter)
    }

    pub fn critical_points_pzetas<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.0.critical_points_pzetas().to_vec()).into_pyarray(py)
    }

    pub fn critical_points_energies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.0.critical_points_energies().to_vec()).into_pyarray(py)
    }

    pub fn boundary_extrema_pzetas<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.0.boundary_extrema_pzetas().to_vec()).into_pyarray(py)
    }

    pub fn boundary_extrema_energies<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        Array1::from_vec(self.0.boundary_extrema_energies().to_vec()).into_pyarray(py)
    }
}

// ===============================================================================================

wrapper_debug_export!(PyParabola);
wrapper_debug_export!(PyEnergyPzetaPlane);
impl_py_repr!(PyParabola, pretty);
impl_py_repr!(PyEnergyPzetaPlane, pretty);
