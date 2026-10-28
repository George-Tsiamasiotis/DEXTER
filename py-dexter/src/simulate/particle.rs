//! Defines the `Particle` wrapper.

use crate::*;
use dexter::dexter_machine::*;
use dexter::dexter_simulate::*;

use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyType;

#[pyclass(name = "_PyParticle", from_py_object)]
#[derive(Clone)]
pub struct PyParticle(pub(crate) Particle);

// ===============================================================================================

#[pymethods] // Routines
impl PyParticle {
    #[new]
    pub fn new(initial: &PyInitialConditions) -> Self {
        Self(Particle::new(&initial.0))
    }

    #[classmethod]
    pub fn from_critical_point(_: &Bound<'_, PyType>, point: &PyCriticalPoint) -> Self {
        Self(Particle::from_critical_point(&point.0))
    }

    pub fn integrate(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        efield: &PyEfield,
        perturbation: &PyPerturbation,
        teval: (f64, f64),
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_efield(efield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.integrate(machine, teval, &solver_params.0);
    }

    pub fn intersect(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        efield: &PyEfield,
        perturbation: &PyPerturbation,
        intersect_params: &PyIntersectParams,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_efield(efield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0
            .intersect(machine, &intersect_params.0, &solver_params.0);
    }

    pub fn close(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        efield: &PyEfield,
        perturbation: &PyPerturbation,
        periods: usize,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_efield(efield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.close(machine, periods, &solver_params.0);
    }

    pub fn classify(&mut self, qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();
        self.0.classify(machine);
    }
}

#[pymethods] // Getters
impl PyParticle {
    #[getter]
    pub fn initial_conditions(&self) -> PyInitialConditions {
        PyInitialConditions(self.0.initial_conditions().clone())
    }

    #[getter]
    pub fn integration_status(&self) -> String {
        format!("{:?}", self.0.integration_status())
    }

    #[getter]
    pub fn steps_taken(&self) -> usize {
        self.0.steps_taken()
    }

    #[getter]
    pub fn steps_stored(&self) -> usize {
        self.0.steps_stored()
    }

    #[getter]
    pub fn duration(&self) -> String {
        format!("{:?}", self.0.duration())
    }

    #[getter]
    pub fn initial_energy(&self) -> Option<f64> {
        self.0.initial_energy()
    }

    #[getter]
    pub fn final_energy(&self) -> Option<f64> {
        self.0.final_energy()
    }

    #[getter]
    pub fn energy_rsd(&self) -> Option<f64> {
        self.0.energy_rsd()
    }

    #[getter]
    pub fn energy_pzeta_position(&self) -> String {
        format!("{:?}", self.0.energy_pzeta_position())
    }

    #[getter]
    pub fn orbit_type(&self) -> String {
        format!("{:?}", self.0.orbit_type())
    }

    #[getter]
    pub fn omega_theta(&self) -> Option<f64> {
        self.0.omega_theta()
    }

    #[getter]
    pub fn omega_zeta(&self) -> Option<f64> {
        self.0.omega_zeta()
    }

    #[getter]
    pub fn qkinetic(&self) -> Option<f64> {
        self.0.qkinetic()
    }

    pub fn print_caches(&self) {
        self.0.print_caches()
    }

    pub fn discard_arrays(&mut self) {
        self.0.discard_arrays()
    }

    #[getter]
    pub fn flux_cache_hits(&self) -> usize {
        dbg!(self.0.flux_cache_hits())
    }

    #[getter]
    pub fn flux_cache_misses(&self) -> usize {
        self.0.flux_cache_misses()
    }

    #[getter]
    pub fn theta_cache_hits(&self) -> usize {
        self.0.theta_cache_hits()
    }

    #[getter]
    pub fn theta_cache_misses(&self) -> usize {
        self.0.theta_cache_misses()
    }

    #[getter]
    pub fn mode_cache_hits(&self) -> usize {
        self.0.mode_cache_hits()
    }

    #[getter]
    pub fn mode_cache_misses(&self) -> usize {
        self.0.mode_cache_misses()
    }

    #[rustfmt::skip]
    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let values: &[f64] = match name {
            "t_array"      => self.0.t_values(),
            "psit_array"   => self.0.psit_values(),
            "psip_array"   => self.0.psip_values(),
            "theta_array"  => self.0.theta_values(),
            "zeta_array"   => self.0.zeta_values(),
            "rho_array"    => self.0.rho_values(),
            "mu_array"     => self.0.mu_values(),
            "ptheta_array" => self.0.ptheta_values(),
            "pzeta_array"  => self.0.pzeta_values(),
            "energy_array" => self.0.energy_values(),
            _ => return Err(DexterError::AttributeError {
                obj: "NcQfactor".into(),
                attr: name.into(),
            }),
        };
        Ok(Array1::from(values.to_vec()).into_pyarray(py))
    }
}

// ===============================================================================================

wrapper_debug_export!(PyParticle);
impl_py_repr!(PyParticle, pretty);
