//! Exports [`EnergyPzetaPlane`] and its helper types.

use dexter::dexter_comspace::*;
use dexter::dexter_machine::*;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use parabola::{Line, LineIntercepts, Parabola};

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

#[pyclass(name = "_PyEnergyPzetaPlane", immutable_type, frozen)]
pub struct PyEnergyPzetaPlane(pub EnergyPzetaPlane);

#[pymethods]
impl PyEnergyPzetaPlane {
    #[new]
    pub fn new(qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield, mu: f64) -> Self {
        // EnergyPzetaPlane does not use a Geometry or a Perturbation
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        Self(EnergyPzetaPlane::from_mu(machine, mu))
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
}

// ===============================================================================================

wrapper_debug_export!(PyParabola);
wrapper_debug_export!(PyEnergyPzetaPlane);
impl_py_repr!(PyParabola, pretty);
impl_py_repr!(PyEnergyPzetaPlane, pretty);
