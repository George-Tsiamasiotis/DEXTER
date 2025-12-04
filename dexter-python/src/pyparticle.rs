//! Particle objects' Python wrappers.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use safe_unwrap::safe_unwrap;

use particle::{Evolution, Frequencies, InitialConditions, MappingParameters, Particle, Radians};
use utils::{
    py_debug_impl, py_export_getter, py_get_enum_string, py_get_numpy1D, py_get_primitive_field,
    py_repr_impl,
};

use crate::pyerrors::PyParticleError;
use crate::pylibrium::{PyBfield, PyCurrents, PyPerturbation, PyQfactor};

#[derive(Clone)]
#[pyclass(frozen, name = "InitialConditions")]
pub struct PyInitialConditions(InitialConditions);

#[pymethods]
impl PyInitialConditions {
    #[new]
    fn new_py(time0: f64, theta0: f64, psip0: f64, rho0: f64, zeta0: f64, mu: f64) -> Self {
        Self(InitialConditions {
            time0,
            theta0,
            psip0,
            rho0,
            zeta0,
            mu,
        })
    }
}

py_debug_impl!(PyInitialConditions);
py_repr_impl!(PyInitialConditions);
py_get_primitive_field!(PyInitialConditions, time0, f64);
py_get_primitive_field!(PyInitialConditions, theta0, f64);
py_get_primitive_field!(PyInitialConditions, psip0, f64);
py_get_primitive_field!(PyInitialConditions, rho0, f64);
py_get_primitive_field!(PyInitialConditions, zeta0, f64);
py_get_primitive_field!(PyInitialConditions, mu, f64);

// ===============================================================================================

#[pyclass(frozen, name = "MappingParameters")]
pub struct PyMappingParameters(pub MappingParameters);

#[pymethods]
impl PyMappingParameters {
    #[new]
    pub fn new_py(section: &str, alpha: Radians, intersections: usize) -> Self {
        let section = match section.to_lowercase().as_str() {
            "theta" | "consttheta" => particle::PoincareSection::ConstTheta,
            "zeta" | "constzeta" => particle::PoincareSection::ConstZeta,
            _ => panic!("Mapping section must be either 'theta' or 'zeta'"),
        };
        Self(MappingParameters::new(section, alpha, intersections))
    }
}

py_debug_impl!(PyMappingParameters);
py_repr_impl!(PyMappingParameters);
py_get_enum_string!(PyMappingParameters, section);
py_get_primitive_field!(PyMappingParameters, alpha, Radians);
py_get_primitive_field!(PyMappingParameters, intersections, usize);

// ===============================================================================================

#[pyclass(frozen, name = "Evolution")]
pub struct PyEvolution(Evolution);

py_debug_impl!(PyEvolution);
py_repr_impl!(PyEvolution);
py_get_numpy1D!(PyEvolution, time);
py_get_numpy1D!(PyEvolution, theta);
py_get_numpy1D!(PyEvolution, psip);
py_get_numpy1D!(PyEvolution, rho);
py_get_numpy1D!(PyEvolution, zeta);
py_get_numpy1D!(PyEvolution, psi);
py_get_numpy1D!(PyEvolution, ptheta);
py_get_numpy1D!(PyEvolution, pzeta);
py_get_numpy1D!(PyEvolution, energy);
py_get_primitive_field!(PyEvolution, energy_std, f64);
py_export_getter!(PyEvolution, steps_taken, usize);
py_export_getter!(PyEvolution, steps_stored, usize);

// ===============================================================================================

#[pyclass(frozen, name = "Frequencies")]
pub struct PyFrequencies(Frequencies);

py_debug_impl!(PyFrequencies);
py_repr_impl!(PyFrequencies);
py_export_getter!(PyFrequencies, omega_theta, Option<f64>);
py_export_getter!(PyFrequencies, omega_zeta, Option<f64>);
py_export_getter!(PyFrequencies, qkinetic, Option<f64>);

// ===============================================================================================

#[pyclass(name = "Particle")]
pub struct PyParticle(pub Particle);

#[pymethods]
impl PyParticle {
    #[new]
    pub fn new_py(initial: PyInitialConditions) -> Self {
        Self(Particle::new(&initial.0))
    }

    #[getter]
    pub fn get_initial_conditions(&self) -> PyInitialConditions {
        PyInitialConditions(self.0.initial_conditions.clone())
    }

    #[getter]
    pub fn get_evolution(&self) -> PyEvolution {
        PyEvolution(self.0.evolution.clone())
    }

    #[getter]
    pub fn get_frequencies(&self) -> PyFrequencies {
        PyFrequencies(self.0.frequencies.clone())
    }

    pub fn integrate<'py>(
        &mut self,
        qfactor: &PyQfactor,
        currents: &PyCurrents,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        t_eval: Bound<'py, PyTuple>,
    ) -> Result<(), PyParticleError> {
        match t_eval.len() {
            2 => (),
            _ => panic!("`t_eval` must be of the form (t0, tf)"),
        };
        let t_eval: Vec<f64> = t_eval
            .iter()
            .map(|any| {
                any.extract::<f64>()
                    .expect("t_eval elements must be floats")
            })
            .collect();
        let t_eval = (
            safe_unwrap!("len already checked", t_eval.first().copied()),
            safe_unwrap!("len already checked", t_eval.last().copied()),
        );

        Ok(self
            .0
            .integrate(&qfactor.0, &bfield.0, &currents.0, &perturbation.0, t_eval)?)
    }

    pub fn map(
        &mut self,
        qfactor: &PyQfactor,
        currents: &PyCurrents,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        params: &PyMappingParameters,
    ) -> Result<(), PyParticleError> {
        Ok(self.0.map(
            &qfactor.0,
            &bfield.0,
            &currents.0,
            &perturbation.0,
            &params.0,
        )?)
    }

    pub fn calculate_frequencies(
        &mut self,
        qfactor: &PyQfactor,
        currents: &PyCurrents,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
    ) -> Result<(), PyParticleError> {
        Ok(self
            .0
            .calculate_frequencies(&qfactor.0, &bfield.0, &currents.0, &perturbation.0)?)
    }
}

py_debug_impl!(PyParticle);
py_repr_impl!(PyParticle);
py_get_enum_string!(PyParticle, status);
py_get_enum_string!(PyParticle, orbit_type);
py_export_getter!(PyParticle, initial_energy, f64);
py_export_getter!(PyParticle, final_energy, f64);
