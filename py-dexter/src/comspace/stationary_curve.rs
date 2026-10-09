use dexter::dexter_comspace::*;
use dexter::dexter_machine::*;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::*;

// ===============================================================================================

#[pyclass(name = "_PyStationaryCurveSegment", immutable_type, frozen)]
pub struct PyStationaryCurveSegment(pub StationaryCurveSegment);

#[pymethods]
impl PyStationaryCurveSegment {
    #[getter]
    pub fn theta<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.theta().to_owned().into_pyarray(py)
    }

    #[getter]
    pub fn flux<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.0.flux().to_owned().into_pyarray(py)
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}

// ===============================================================================================

#[pyclass(name = "_PyStationaryCurve", immutable_type, frozen)]
pub struct PyStationaryCurve(pub StationaryCurve);

#[pymethods]
impl PyStationaryCurve {
    #[new]
    pub fn new(qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield) -> Result<Self> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        Ok(Self(StationaryCurve::build(machine)?))
    }

    #[getter]
    pub fn flux_kind(&self) -> String {
        format!("{:?}", self.0.flux_kind)
    }

    #[getter]
    pub fn segments<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyList>> {
        let py_segments: Vec<PyStationaryCurveSegment> = self
            .0
            .segments
            .iter()
            .map(|segment| PyStationaryCurveSegment(segment.clone()))
            .collect();
        Ok(PyList::new(py, py_segments)?)
    }
}

// ===============================================================================================

wrapper_debug_export!(PyStationaryCurveSegment);
wrapper_debug_export!(PyStationaryCurve);
impl_py_repr!(PyStationaryCurveSegment, pretty);
impl_py_repr!(PyStationaryCurve, pretty);
