//! `dexter_simulate::Queue` newtype, constructors and method exports.

use dexter::dexter_simulate::{
    FluxCoordinate, Queue, QueueInitialConditions, poloidal_fluxes, toroidal_fluxes,
};
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::pyerror::PySimulationError;
use crate::{
    generic_queue_classify_impl, py_debug_impl, py_get_enum_string, py_get_numpy1D,
    py_get_numpy1D_fallible, py_repr_impl,
};

use crate::pyparticle::{PyInitialConditions, PyParticle};

// ===============================================================================================

#[pyclass(name = "_PyInitialFluxArray", frozen)]
pub struct PyInitialFluxArray {
    kind: FluxCoordinate,
    values: Array1<f64>,
}

#[pymethods]
impl PyInitialFluxArray {
    #[new]
    #[pyo3(signature = (kind, values))]
    pub fn new(kind: &str, values: Vec<f64>) -> PyResult<Self> {
        let kind = match kind.to_lowercase().as_str() {
            "toroidal" => FluxCoordinate::Toroidal,
            "poloidal" => FluxCoordinate::Poloidal,
            _ => return Err(PyErr::new::<PyTypeError, _>("Invalid 'kind'")),
        };
        let values = Array1::from_vec(values);
        Ok(Self { kind, values })
    }

    #[getter]
    pub fn get_kind(&self) -> String {
        format!("{:?}", self.kind)
    }

    #[getter]
    pub fn get_values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.values.clone().into_pyarray(py)
    }
}

// ===============================================================================================

#[pyclass(name = "_PyQueueInitialConditions", frozen)]
pub struct PyQueueInitialConditions(pub(crate) QueueInitialConditions);

#[pymethods]
impl PyQueueInitialConditions {
    #[classmethod]
    #[pyo3(signature = (t0, flux0, theta0, zeta0, rho0, mu0))]
    pub fn boozer(
        _: &Bound<'_, PyType>,
        t0: Vec<f64>,
        flux0: &PyInitialFluxArray,
        theta0: Vec<f64>,
        zeta0: Vec<f64>,
        rho0: Vec<f64>,
        mu0: Vec<f64>,
    ) -> Result<Self, PySimulationError> {
        let flux0_values = flux0
            .values
            .as_slice()
            .expect("should be in standard layout");
        let flux0_rust = match flux0.kind {
            FluxCoordinate::Toroidal => toroidal_fluxes(flux0_values).to_vec(),
            FluxCoordinate::Poloidal => poloidal_fluxes(flux0_values).to_vec(),
        };
        Ok(Self(QueueInitialConditions::boozer(
            &t0,
            &flux0_rust,
            &theta0,
            &zeta0,
            &rho0,
            &mu0,
        )?))
    }

    #[classmethod]
    #[pyo3(signature = (t0, flux0, theta0, zeta0, pzeta0, mu0))]
    pub fn mixed(
        _: &Bound<'_, PyType>,
        t0: Vec<f64>,
        flux0: &PyInitialFluxArray,
        theta0: Vec<f64>,
        zeta0: Vec<f64>,
        pzeta0: Vec<f64>,
        mu0: Vec<f64>,
    ) -> Result<Self, PySimulationError> {
        let flux0_values = flux0
            .values
            .as_slice()
            .expect("should be in standard layout");
        let flux0_rust = match flux0.kind {
            FluxCoordinate::Toroidal => toroidal_fluxes(flux0_values).to_vec(),
            FluxCoordinate::Poloidal => poloidal_fluxes(flux0_values).to_vec(),
        };
        Ok(Self(QueueInitialConditions::mixed(
            &t0,
            &flux0_rust,
            &theta0,
            &zeta0,
            &pzeta0,
            &mu0,
        )?))
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }

    pub fn __getitem__(&self, index: usize) -> PyInitialConditions {
        PyInitialConditions(self.0.initial_from_index(index).clone())
    }
}

py_debug_impl!(PyQueueInitialConditions);
py_repr_impl!(PyQueueInitialConditions);
py_get_numpy1D!(PyQueueInitialConditions, t_array);
py_get_numpy1D!(PyQueueInitialConditions, theta_array);
py_get_numpy1D!(PyQueueInitialConditions, zeta_array);
py_get_numpy1D!(PyQueueInitialConditions, mu_array);
py_get_numpy1D_fallible!(PyQueueInitialConditions, rho_array);
py_get_numpy1D_fallible!(PyQueueInitialConditions, pzeta_array);

// ===============================================================================================

#[pyclass(name = "_PyQueue")]
pub struct PyQueue(pub(crate) Queue);

#[pymethods]
impl PyQueue {
    #[new]
    #[pyo3(signature = (initial_conditions))]
    pub fn new(initial_conditions: &PyQueueInitialConditions) -> Self {
        Self(Queue::new(&initial_conditions.0))
    }

    #[getter]
    pub fn get_initial_conditions(&self) -> PyQueueInitialConditions {
        PyQueueInitialConditions(self.0.initial_conditions())
    }

    #[getter]
    pub fn get_particle_count(&self) -> usize {
        self.0.particle_count()
    }

    #[getter]
    /// Creates a python list with the contained particles
    pub fn get_particles(&self) -> Vec<PyParticle> {
        self.0
            .iter()
            .map(|particle| PyParticle(particle.clone()))
            .collect()
    }

    #[getter]
    pub fn get_steps_taken_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.0.steps_taken_array().into_pyarray(py)
    }

    #[getter]
    pub fn get_steps_stored_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.0.steps_stored_array().into_pyarray(py)
    }

    #[getter("_durations_as_nanos")]
    pub fn get_durations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        Array1::from_iter(
            self.0
                .durations()
                .iter()
                .map(|duration| duration.as_nanos() as usize),
        )
        .into_pyarray(py)
    }

    pub fn __getitem__(&self, index: usize) -> PyParticle {
        PyParticle(self.0[index].clone())
    }
}

py_debug_impl!(PyQueue);
py_repr_impl!(PyQueue);
py_get_enum_string!(PyQueue, routine);
py_get_numpy1D!(PyQueue, energy_array);
py_get_numpy1D!(PyQueue, omega_theta_array);
py_get_numpy1D!(PyQueue, omega_zeta_array);
py_get_numpy1D!(PyQueue, qkinetic_array);

#[allow(non_camel_case_types)]
/// Since we can not use generics in top level types and methods, we must do the compiler's work
/// and do a manual monomorphization.
///
/// Since #[pymethods] restricts the use of macros, this is the best I could come up with.
///
/// The naming format is
///
///     `__<method>_<Qfactor._dyn>_<Current._dyn>_<Bfield._dyn>_<Perturbation._dyn>`
///
/// where `_dyn` is a 4 letter class attribute string in each final wrapped object. Each
/// routine's wrapper concatenates these attributes to form the monomorphized methods name and
/// calls it.
mod py_particle_generics_impl {
    use super::*;
    use dexter::dexter_simulate::SolverParams;

    use crate::pyparticle::{PyIntersectParams, resolve_stepping_method};

    use crate::{
        generic_queue_close_impl, generic_queue_integrate_impl, generic_queue_intersect_impl,
    };

    type uniQ = crate::pyqfactors::PyUnityQfactor;
    type parQ = crate::pyqfactors::PyParabolicQfactor;
    type ncdQ = crate::pyqfactors::PyNcQfactor;

    type larC = crate::pycurrents::PyLarCurrent;
    type ncdC = crate::pycurrents::PyNcCurrent;

    type larB = crate::pybfields::PyLarBfield;
    type ncdB = crate::pybfields::PyNcBfield;

    type cosP = crate::pyperturbation::PyCosPerturbation;
    type ncdP = crate::pyperturbation::PyNcPerturbation;

    // ================= Integrate

    generic_queue_integrate_impl!(__integrate_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_queue_integrate_impl!(__integrate_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_queue_integrate_impl!(__integrate_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_queue_integrate_impl!(__integrate_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_queue_integrate_impl!(__integrate_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_queue_integrate_impl!(__integrate_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_queue_integrate_impl!(__integrate_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_queue_integrate_impl!(__integrate_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_queue_integrate_impl!(__integrate_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);

    // ================= Intersect

    generic_queue_intersect_impl!(__intersect_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_queue_intersect_impl!(__intersect_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_queue_intersect_impl!(__intersect_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_queue_intersect_impl!(__intersect_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_queue_intersect_impl!(__intersect_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_queue_intersect_impl!(__intersect_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_queue_intersect_impl!(__intersect_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_queue_intersect_impl!(__intersect_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_queue_intersect_impl!(__intersect_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);

    // ================= Close

    generic_queue_close_impl!(__close_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_queue_close_impl!(__close_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_queue_close_impl!(__close_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_queue_close_impl!(__close_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_queue_close_impl!(__close_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_queue_close_impl!(__close_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_queue_close_impl!(__close_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_queue_close_impl!(__close_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_queue_close_impl!(__close_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_queue_close_impl!(__close_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_queue_close_impl!(__close_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_queue_close_impl!(__close_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_queue_close_impl!(__close_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_queue_close_impl!(__close_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_queue_close_impl!(__close_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_queue_close_impl!(__close_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_queue_close_impl!(__close_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_queue_close_impl!(__close_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_queue_close_impl!(__close_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_queue_close_impl!(__close_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_queue_close_impl!(__close_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_queue_close_impl!(__close_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_queue_close_impl!(__close_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_queue_close_impl!(__close_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);

    // ================= Classify

    generic_queue_classify_impl!(__classify_uniQ_larC_larB, uniQ, larC, larB);
    generic_queue_classify_impl!(__classify_uniQ_larC_ncdB, uniQ, larC, ncdB);
    generic_queue_classify_impl!(__classify_uniQ_ncdC_larB, uniQ, ncdC, larB);
    generic_queue_classify_impl!(__classify_uniQ_ncdC_ncdB, uniQ, ncdC, ncdB);

    generic_queue_classify_impl!(__classify_parQ_larC_larB, parQ, larC, larB);
    generic_queue_classify_impl!(__classify_parQ_larC_ncdB, parQ, larC, ncdB);
    generic_queue_classify_impl!(__classify_parQ_ncdC_larB, parQ, ncdC, larB);
    generic_queue_classify_impl!(__classify_parQ_ncdC_ncdB, parQ, ncdC, ncdB);

    generic_queue_classify_impl!(__classify_ncdQ_larC_larB, ncdQ, larC, larB);
    generic_queue_classify_impl!(__classify_ncdQ_larC_ncdB, ncdQ, larC, ncdB);
    generic_queue_classify_impl!(__classify_ncdQ_ncdC_larB, ncdQ, ncdC, larB);
    generic_queue_classify_impl!(__classify_ncdQ_ncdC_ncdB, ncdQ, ncdC, ncdB);
}
