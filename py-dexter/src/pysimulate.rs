//! `dexter-simulate` newtypes, constructors and method exports.

use dexter::dexter_simulate::{
    FluxCoordinate, InitialConditions, InitialFlux, IntersectParams, Intersection, Particle, Queue,
    QueueInitialConditions, SolverParams, SteppingMethod,
};
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::pyerror::PySimulationError;
use crate::pylibrium::*;
use crate::{
    generic_particle_integrate_impl, generic_particle_intersect_impl, generic_queue_integrate_impl,
    generic_queue_intersect_impl, py_debug_impl, py_export_getter, py_export_pub_field,
    py_get_enum_string, py_get_numpy1D, py_repr_impl,
};

// ===============================================================================================

#[pyclass(name = "_PyInitialFlux", frozen)]
pub struct PyInitialFlux(pub(crate) InitialFlux);

#[pymethods]
impl PyInitialFlux {
    #[new]
    #[pyo3(signature = (kind, value))]
    pub fn new(kind: &str, value: f64) -> PyResult<Self> {
        match kind.to_lowercase().as_str() {
            "toroidal" => Ok(Self(InitialFlux::Toroidal(value))),
            "poloidal" => Ok(Self(InitialFlux::Poloidal(value))),
            _ => Err(PyErr::new::<PyTypeError, _>("Invalid 'InitialFlux'")),
        }
    }

    #[getter]
    pub fn get_kind(&self) -> &str {
        match self.0 {
            InitialFlux::Toroidal(_) => "Toroidal",
            InitialFlux::Poloidal(_) => "Poloidal",
        }
    }

    #[getter]
    pub fn get_value(&self) -> f64 {
        self.0.value()
    }
}

py_debug_impl!(PyInitialFlux);
py_repr_impl!(PyInitialFlux);

// ===============================================================================================
#[derive(Clone)]
#[pyclass(name = "_PyInitialConditions", frozen)]
pub struct PyInitialConditions(pub(crate) InitialConditions);

#[pymethods]
impl PyInitialConditions {
    #[new]
    pub fn new<'py>(
        t0: f64,
        flux0: &PyInitialFlux,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> PyResult<Self> {
        let flux0: InitialFlux = flux0.0;
        Ok(Self(InitialConditions::boozer(
            t0, flux0, theta0, zeta0, rho0, mu0,
        )))
    }

    #[getter]
    pub fn get_flux0(&self) -> PyInitialFlux {
        PyInitialFlux(self.0.flux0)
    }
}

py_debug_impl!(PyInitialConditions);
py_repr_impl!(PyInitialConditions);
py_export_pub_field!(PyInitialConditions, t0, f64);
py_export_pub_field!(PyInitialConditions, theta0, f64);
py_export_pub_field!(PyInitialConditions, zeta0, f64);
py_export_pub_field!(PyInitialConditions, rho0, f64);
py_export_pub_field!(PyInitialConditions, mu0, f64);

// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyIntersectParams", frozen)]
pub struct PyIntersectParams(pub(crate) IntersectParams);

#[pymethods]
impl PyIntersectParams {
    #[new]
    pub fn new<'py>(intersection: &str, angle: f64, turns: usize) -> PyResult<Self> {
        let intersection = match intersection.to_lowercase().as_str() {
            "consttheta" => Intersection::ConstTheta,
            "constzeta" => Intersection::ConstZeta,
            _ => {
                return Err(PyErr::new::<PyTypeError, _>(
                    "'intersection' must be either 'ConstTheta' or 'ConstZeta'",
                ));
            }
        };
        Ok(Self(IntersectParams::new(intersection, angle, turns)))
    }

    #[getter]
    pub fn intersection(&self) -> String {
        format!("{:#?}", self.0.intersection)
    }
}

py_debug_impl!(PyIntersectParams);
py_repr_impl!(PyIntersectParams);
py_export_pub_field!(PyIntersectParams, angle, f64);
py_export_pub_field!(PyIntersectParams, turns, usize);

// ===============================================================================================

#[pyclass(name = "_PyParticle")]
pub struct PyParticle(pub(crate) Particle);

#[pymethods]
impl PyParticle {
    #[new]
    #[pyo3(signature = (initial_conditions))]
    pub fn new(initial_conditions: PyInitialConditions) -> Self {
        Self(Particle::new(&initial_conditions.0))
    }

    #[getter]
    pub fn get_initial_conditions(&self) -> PyInitialConditions {
        PyInitialConditions(self.0.initial_conditions())
    }

    pub fn print_stats(&self) {
        self.0.print_cache_stats();
    }
}

/// Attempts to create a valid [`SteppingMethod`] object from the corresponding [`PyAny`] Python
/// argument.
///
/// Returns early if a valid string is found that matches one of the simple variants,
/// otherwise tries to unpack the [`PyAny`] object into a ("FixedStep", f64) tuple and cast it
/// into a SteppingMethod::FixedStep.
///
/// Returns an Error if both string matching and casting fail.
fn resolve_stepping_method<'py>(arg: Bound<'py, PyAny>) -> PyResult<SteppingMethod> {
    use SteppingMethod::*;
    match arg.to_string().to_lowercase().as_str() {
        "energyadaptivestep" => return Ok(EnergyAdaptiveStep),
        "erroradaptivestep" => return Ok(ErrorAdaptiveStep),
        _ => (),
    };
    let tuple = arg.cast::<PyTuple>()?;
    let string: String = tuple.get_item(0)?.extract::<String>()?.to_lowercase();
    let value: f64 = tuple.get_item(1)?.extract::<f64>()?;
    match string.as_str() {
        "fixedstep" if value.is_finite() => Ok(FixedStep(value)),
        _ => Err(PyErr::new::<PyTypeError, _>("Invalid 'method'")),
    }
}

py_debug_impl!(PyParticle);
py_repr_impl!(PyParticle);
py_get_enum_string!(PyParticle, integration_status);
py_export_getter!(PyParticle, steps_taken, usize);
py_export_getter!(PyParticle, steps_stored, usize);
py_export_getter!(PyParticle, initial_energy, Option<f64>);
py_export_getter!(PyParticle, final_energy, Option<f64>);
py_export_getter!(PyParticle, energy_var, Option<f64>);
py_get_numpy1D!(PyParticle, t_array);
py_get_numpy1D!(PyParticle, psi_array);
py_get_numpy1D!(PyParticle, psip_array);
py_get_numpy1D!(PyParticle, theta_array);
py_get_numpy1D!(PyParticle, zeta_array);
py_get_numpy1D!(PyParticle, rho_array);
py_get_numpy1D!(PyParticle, mu_array);
py_get_numpy1D!(PyParticle, ptheta_array);
py_get_numpy1D!(PyParticle, pzeta_array);
py_get_numpy1D!(PyParticle, energy_array);

// ===============================================================================================

#[pyclass(name = "_PyInitialFluxArray1", frozen)]
pub struct PyInitialFluxArray1 {
    kind: FluxCoordinate,
    values: Array1<f64>,
}

#[pymethods]
impl PyInitialFluxArray1 {
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
    #[new]
    #[pyo3(signature = (t0, flux0, theta0, zeta0, rho0, mu0))]
    pub fn new<'py>(
        t0: Vec<f64>,
        flux0: &PyInitialFluxArray1,
        theta0: Vec<f64>,
        zeta0: Vec<f64>,
        rho0: Vec<f64>,
        mu0: Vec<f64>,
    ) -> Result<Self, PySimulationError> {
        let flux0_rust: Array1<InitialFlux>;
        match flux0.kind {
            FluxCoordinate::Toroidal => {
                flux0_rust =
                    Array1::from_iter(flux0.values.map(|value| InitialFlux::Toroidal(*value)))
            }
            FluxCoordinate::Poloidal => {
                flux0_rust =
                    Array1::from_iter(flux0.values.map(|value| InitialFlux::Poloidal(*value)))
            }
        };
        Ok(Self(QueueInitialConditions::build(
            &t0,
            &flux0_rust.as_slice().unwrap_or_default(), // let rust handle it
            &theta0,
            &zeta0,
            &rho0,
            &mu0,
        )?))
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}

py_debug_impl!(PyQueueInitialConditions);
py_repr_impl!(PyQueueInitialConditions);
py_get_numpy1D!(PyQueueInitialConditions, t_array);
py_get_numpy1D!(PyQueueInitialConditions, theta_array);
py_get_numpy1D!(PyQueueInitialConditions, zeta_array);
py_get_numpy1D!(PyQueueInitialConditions, rho_array);
py_get_numpy1D!(PyQueueInitialConditions, mu_array);

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
    pub fn initial_conditions(&self) -> PyQueueInitialConditions {
        PyQueueInitialConditions(self.0.initial_conditions())
    }

    #[getter]
    pub fn particle_count(&self) -> usize {
        self.0.particle_count()
    }

    #[getter]
    /// Creates a python list with the contained particles
    pub fn particles(&self) -> Vec<PyParticle> {
        self.0
            .iter()
            .map(|particle| PyParticle(particle.clone()))
            .collect()
    }

    pub fn __getitem__(&self, index: usize) -> PyParticle {
        PyParticle(self.0[index].clone())
    }
}

py_debug_impl!(PyQueue);
py_repr_impl!(PyQueue);
py_get_enum_string!(PyQueue, routine);

// ===============================================================================================
// ===============================================================================================
// ===============================================================================================

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

    type uniQ = PyUnityQfactor;
    type parQ = PyParabolicQfactor;
    type ncdQ = PyNcQfactor;

    type larC = PyLarCurrent;
    type ncdC = PyNcCurrent;

    type larB = PyLarBfield;
    type ncdB = PyNcBfield;

    type cosP = PyCosPerturbation;
    type ncdP = PyNcPerturbation;

    // ================= Integrate

    generic_particle_integrate_impl!(__integrate_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__integrate_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_particle_integrate_impl!(__integrate_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__integrate_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_particle_integrate_impl!(__integrate_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__integrate_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__integrate_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__integrate_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__integrate_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);

    // ================= Intersect

    generic_particle_intersect_impl!(__intersect_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_particle_intersect_impl!(__intersect_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_particle_intersect_impl!(__intersect_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_particle_intersect_impl!(__intersect_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_particle_intersect_impl!(__intersect_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_particle_intersect_impl!(__intersect_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_particle_intersect_impl!(__intersect_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_particle_intersect_impl!(__intersect_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_particle_intersect_impl!(__intersect_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);

    // ================= Integrate Queue

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

    // ================= Intersect Queue

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
}
