//! Defines the `GCHamiltonian` wrapper.

use crate::*;
use dexter::dexter_comspace::*;
use dexter::dexter_machine::*;

use numpy::PyReadonlyArray1;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyType};

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
    pub fn eval_array<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_array: PyReadonlyArray1<f64>,
        flux_array: &PyMagneticFluxArray,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray1<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let array =
            GCHamiltonian::eval_array(machine, &theta_array.as_array(), &flux_array.0, pzeta, mu)?;
        Ok(array.into_pyarray(py))
    }

    #[classmethod]
    pub fn eval_grid<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_array: PyReadonlyArray1<f64>,
        flux_array: &PyMagneticFluxArray,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let grid =
            GCHamiltonian::eval_grid(machine, &theta_array.as_array(), &flux_array.0, pzeta, mu)?;
        Ok(grid.into_pyarray(py))
    }

    #[classmethod]
    pub fn eval_theta_dot_array<'py>(
        _: &Bound<'_, PyType>,
        py: Python<'py>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        theta_array: PyReadonlyArray1<f64>,
        flux_array: &PyMagneticFluxArray,
        pzeta: f64,
        mu: f64,
    ) -> Result<Bound<'py, PyArray1<f64>>> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let array = GCHamiltonian::eval_theta_dot_array(
            machine,
            &theta_array.as_array(),
            &flux_array.0,
            pzeta,
            mu,
        )?;
        Ok(array.into_pyarray(py))
    }
}

// ===============================================================================================

#[pyclass(name = "_PyStationaryCurve", immutable_type, frozen)]
pub struct PyStationaryCurve(pub StationaryCurve);

#[pymethods]
impl PyStationaryCurve {
    #[new]
    pub fn new(qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield) -> Self {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        Self(StationaryCurve::new(machine))
    }

    #[getter]
    pub fn segments<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyList>> {
        let py_segments: Vec<PyStationaryCurveSegment> = self
            .0
            .segments()
            .iter()
            .map(|segment| PyStationaryCurveSegment(segment.clone()))
            .collect();
        Ok(PyList::new(py, py_segments)?)
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}

// ===============================================================================================

#[pyclass(name = "_PyStationaryCurveSegment", immutable_type, frozen)]
pub struct PyStationaryCurveSegment(pub StationaryCurveSegment);

#[pymethods]
impl PyStationaryCurveSegment {
    #[getter]
    pub fn flux_kind(&self) -> String {
        format!("{:?}", self.0.flux_array().kind())
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }

    pub fn flux_array(&self) -> PyMagneticFluxArray {
        PyMagneticFluxArray(self.0.flux_array().to_owned())
    }

    pub fn theta_array<'py>(&self, py: Python<'py>) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        Ok(Some(self.0.theta_array().to_owned().into_pyarray(py)))
    }
}

// ===============================================================================================

#[pyclass(name = "_PyCriticalPoint", immutable_type, frozen)]
pub struct PyCriticalPoint(pub CriticalPoint);

#[pymethods]
impl PyCriticalPoint {
    #[getter]
    pub fn theta(&self) -> f64 {
        self.0.theta()
    }

    #[getter]
    pub fn flux(&self) -> PyMagneticFlux {
        PyMagneticFlux(self.0.flux())
    }

    #[getter]
    pub fn energy(&self) -> f64 {
        self.0.energy()
    }

    #[getter]
    pub fn pzeta(&self) -> f64 {
        self.0.pzeta()
    }

    #[getter]
    pub fn mu(&self) -> f64 {
        self.0.mu()
    }

    #[getter]
    pub fn kind(&self) -> String {
        format!("{:?}", self.0.kind())
    }
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

wrapper_debug_export!(PyStationaryCurve);
impl_py_repr!(PyStationaryCurve, pretty);
wrapper_debug_export!(PyStationaryCurveSegment);
impl_py_repr!(PyStationaryCurveSegment, pretty);
wrapper_debug_export!(PyCriticalPoint);
impl_py_repr!(PyCriticalPoint, pretty);
