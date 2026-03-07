//! `dexter-simulate` newtypes, constructors and method exports.

use dexter::dexter_simulate::{
    InitialConditions, InitialFlux, IntersectParams, Intersection, Particle, SolverParams,
    SteppingMethod,
};
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::pylibrium::*;
use crate::{
    generic_particle_integrate_impl, py_debug_impl, py_export_getter, py_export_pub_field,
    py_get_enum_string, py_get_numpy1D, py_repr_impl,
};

// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyInitialConditions", subclass, frozen)]
pub struct PyInitialConditions(pub(crate) InitialConditions);

#[pymethods]
impl PyInitialConditions {
    #[new]
    pub fn new<'py>(
        t0: f64,
        flux0: Bound<'py, PyTuple>,
        theta0: f64,
        zeta0: f64,
        rho0: f64,
        mu0: f64,
    ) -> PyResult<Self> {
        let tuple = flux0.cast::<PyTuple>()?;
        let flux: String = tuple.get_item(0)?.extract::<String>()?.to_lowercase();
        let value: f64 = tuple.get_item(1)?.extract::<f64>()?;
        let flux0 = match flux.as_str() {
            "toroidal" if value.is_finite() => InitialFlux::Toroidal(value),
            "poloidal" if value.is_finite() => InitialFlux::Poloidal(value),
            _ => return Err(PyErr::new::<PyTypeError, _>("Invalid 'flux0'")),
        };
        Ok(Self(InitialConditions {
            t0,
            flux0,
            theta0,
            zeta0,
            rho0,
            mu0,
        }))
    }

    #[getter]
    pub fn get_flux0<'py>(&self, py: Python<'py>) -> Result<Bound<'py, PyTuple>, PyErr> {
        use InitialFlux::*;
        let tuple: (String, f64);
        match self.0.flux0 {
            Toroidal(psi) => tuple = ("Toroidal".into(), psi),
            Poloidal(psip) => tuple = ("Poloidal".into(), psip),
        }
        tuple.into_pyobject(py)
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
#[pyclass(name = "_PyIntersectParams", subclass, frozen)]
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

#[pyclass(name = "_PyParticle", subclass)]
pub struct PyParticle(pub(crate) Particle);

#[pymethods]
impl PyParticle {
    #[new]
    pub fn new(initial: PyInitialConditions) -> Self {
        Self(Particle::new(&initial.0))
    }

    #[getter]
    pub fn initial_conditions(&self) -> PyInitialConditions {
        PyInitialConditions(self.0.initial_conditions().clone())
    }
}

impl PyParticle {
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

    generic_particle_integrate_impl!(__int_uniQ_larC_larB_cosP, uniQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__int_uniQ_larC_larB_ncdP, uniQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__int_uniQ_larC_ncdB_cosP, uniQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_uniQ_larC_ncdB_ncdP, uniQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__int_uniQ_ncdC_larB_cosP, uniQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__int_uniQ_ncdC_larB_ncdP, uniQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__int_uniQ_ncdC_ncdB_cosP, uniQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_uniQ_ncdC_ncdB_ncdP, uniQ, ncdC, ncdB, ncdP);

    generic_particle_integrate_impl!(__int_parQ_larC_larB_cosP, parQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__int_parQ_larC_larB_ncdP, parQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__int_parQ_larC_ncdB_cosP, parQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_parQ_larC_ncdB_ncdP, parQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__int_parQ_ncdC_larB_cosP, parQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__int_parQ_ncdC_larB_ncdP, parQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__int_parQ_ncdC_ncdB_cosP, parQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_parQ_ncdC_ncdB_ncdP, parQ, ncdC, ncdB, ncdP);

    generic_particle_integrate_impl!(__int_ncdQ_larC_larB_cosP, ncdQ, larC, larB, cosP);
    generic_particle_integrate_impl!(__int_ncdQ_larC_larB_ncdP, ncdQ, larC, larB, ncdP);
    generic_particle_integrate_impl!(__int_ncdQ_larC_ncdB_cosP, ncdQ, larC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_ncdQ_larC_ncdB_ncdP, ncdQ, larC, ncdB, ncdP);
    generic_particle_integrate_impl!(__int_ncdQ_ncdC_larB_cosP, ncdQ, ncdC, larB, cosP);
    generic_particle_integrate_impl!(__int_ncdQ_ncdC_larB_ncdP, ncdQ, ncdC, larB, ncdP);
    generic_particle_integrate_impl!(__int_ncdQ_ncdC_ncdB_cosP, ncdQ, ncdC, ncdB, cosP);
    generic_particle_integrate_impl!(__int_ncdQ_ncdC_ncdB_ncdP, ncdQ, ncdC, ncdB, ncdP);
}
