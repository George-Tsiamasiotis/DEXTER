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

    #[classmethod]
    pub fn build_stationary_curve(
        _: &Bound<'_, PyType>,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
    ) -> PyStationaryCurve {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        let curve = GCHamiltonian::build_stationary_curve(machine);
        PyStationaryCurve(curve)
    }
}

// ===============================================================================================

#[pyclass(name = "_PyStationaryCurve", immutable_type, frozen)]
pub struct PyStationaryCurve(pub StationaryCurve);

#[pymethods]
impl PyStationaryCurve {
    #[getter]
    pub fn segments<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyList>> {
        let py_segments: Vec<PyParametrizedSegment> = self
            .0
            .segments()
            .iter()
            .map(|segment| PyParametrizedSegment(segment.clone()))
            .collect();
        Ok(PyList::new(py, py_segments)?)
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}

// ===============================================================================================

#[pyclass(name = "_PyParametrizedSegment")]
pub struct PyParametrizedSegment(pub ParametrizedSegment);

#[pymethods]
impl PyParametrizedSegment {
    #[getter]
    pub fn flux_kind(&self) -> String {
        format!("{:?}", self.0.flux_array().kind())
    }

    pub fn calculate_theta_dot(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        pzeta: f64,
        mu: f64,
    ) -> Result<()> {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        self.0.calculate_theta_dot(machine, pzeta, mu)?;
        Ok(())
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }

    pub fn flux_array(&self) -> PyMagneticFluxArray {
        PyMagneticFluxArray(self.0.flux_array().to_owned())
    }

    pub fn get_array<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        let array = match name {
            "t_array" => self.0.t_array(),
            "theta_array" => self.0.theta_array(),
            "theta_dot_array" => match self.0.theta_dot_array() {
                Some(array) => array,
                None => return Ok(None),
            },
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "ParametrizedSegment".into(),
                    attr: name.into(),
                });
            }
        };
        Ok(Some(array.to_owned().into_pyarray(py)))
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
wrapper_debug_export!(PyParametrizedSegment);
impl_py_repr!(PyParametrizedSegment, pretty);
