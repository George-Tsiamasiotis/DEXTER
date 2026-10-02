//! Defines the wrapper around [`dexter::dexter_simulate::Queue`].

use dexter::dexter_machine::*;
use dexter::dexter_simulate::*;
use numpy::{IntoPyArray, PyArray1};
use pyo3::types::{PyList, PyType};

use crate::*;
use pyo3::prelude::*;

// ===============================================================================================

#[pyclass(name = "_PyQueue", immutable_type)]
pub struct PyQueue(Queue);

#[pymethods]
impl PyQueue {
    #[new]
    pub fn new<'py>(initial: &PyQueueInitialConditions) -> Self {
        Self(Queue::new(&initial.0))
    }

    #[classmethod]
    pub fn from_particles<'py>(
        _: &Bound<'_, PyType>,
        pyparticles: Bound<'py, PyList>,
    ) -> PyResult<Self> {
        let mut particles = Vec::<Particle>::with_capacity(pyparticles.len());
        for pyparticle in pyparticles {
            particles.push(pyparticle.extract::<PyParticle>()?.0);
        }
        Ok(Self(Queue::from_particles(&particles)))
    }
}

#[pymethods] // Routines
impl PyQueue {
    pub fn integrate(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        teval: (f64, f64),
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.integrate(machine, teval, &solver_params.0);
    }

    pub fn intersect(
        &mut self,
        qfactor: &PyQfactor,
        current: &PyCurrent,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        intersect_params: &PyIntersectParams,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
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
        perturbation: &PyPerturbation,
        periods: usize,
        discard_arrays: bool,
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0
            .close(machine, periods, discard_arrays, &solver_params.0);
    }

    pub fn classify(&mut self, qfactor: &PyQfactor, current: &PyCurrent, bfield: &PyBfield) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner()).build();

        let mu0 = self
            .0
            .initial_conditions()
            .mu_array()
            .first()
            .copied()
            .expect("mu array is non-empty by construction");
        if self
            .0
            .initial_conditions()
            .mu_array()
            .iter()
            .all(|mu| mu == &mu0)
        {
            self.0.classify_common_mu(machine);
        } else {
            self.0.classify(machine);
        }
    }
}

#[pymethods] // Slicing
impl PyQueue {
    pub fn retain_pzeta<'py>(&mut self, start: f64, end: f64) {
        self.0.retain_pzeta(start..end);
    }

    pub fn retain_energy<'py>(&mut self, start: f64, end: f64) {
        self.0.retain_energy(start..end);
    }

    pub fn retain_energy_pzeta_positions(&mut self, positions: Vec<String>) {
        let positions: Vec<EnergyPzetaPosition> = positions
            .iter()
            .map(|pos| match pos.as_str() {
                "Alpha" => EnergyPzetaPosition::Alpha,
                "Beta" => EnergyPzetaPosition::Beta,
                "Gamma" => EnergyPzetaPosition::Gamma,
                "Delta" => EnergyPzetaPosition::Delta,
                "Epsilon" => EnergyPzetaPosition::Epsilon,
                "Zeta" => EnergyPzetaPosition::Zeta,
                "Eta" => EnergyPzetaPosition::Eta,
                "Theta" => EnergyPzetaPosition::Theta,
                "Iota" => EnergyPzetaPosition::Iota,
                "Kappa" => EnergyPzetaPosition::Kappa,
                "Lambda" => EnergyPzetaPosition::Lambda,
                "Mu" => EnergyPzetaPosition::Mu,
                "Undefined" => EnergyPzetaPosition::Undefined,
                "Forbidden" => EnergyPzetaPosition::Forbidden,
                _ => panic!("'positions' must be valid 'EnergyPzetaPosition' variants"),
            })
            .collect();

        self.0.retain_energy_pzeta_positions(&positions);
    }

    pub fn retain_orbit_types(&mut self, orbit_types: Vec<String>) {
        let orbit_types: Vec<OrbitType> = orbit_types
            .iter()
            .map(|typ| match typ.as_str() {
                "Unclassified" => OrbitType::Unclassified,
                "TrappedLost" => OrbitType::TrappedLost,
                "TrappedConfined" => OrbitType::TrappedConfined,
                "CoPassingLost" => OrbitType::CoPassingLost,
                "CoPassingConfined" => OrbitType::CoPassingConfined,
                "CuPassingLost" => OrbitType::CuPassingLost,
                "CuPassingConfined" => OrbitType::CuPassingConfined,
                "Potato" => OrbitType::Potato,
                "Stagnated" => OrbitType::Stagnated,
                "Undefined" => OrbitType::Undefined,
                _ => panic!("'orbit_types' must be valid 'OrbitType' variants"),
            })
            .collect();

        self.0.retain_orbit_types(&orbit_types);
    }
}

#[pymethods] // Getters
impl PyQueue {
    #[getter]
    pub fn initial_conditions(&mut self) -> PyQueueInitialConditions {
        PyQueueInitialConditions(self.0.initial_conditions().clone())
    }

    #[getter]
    pub fn routines<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let routine_strings: Vec<String> = self
            .0
            .routines()
            .iter()
            .map(|routine| format!("{:?}", routine))
            .collect();
        PyList::new(py, routine_strings.iter())
    }

    pub fn particles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.0
                .particles()
                .iter()
                .map(|particle| PyParticle(particle.clone()))
                .collect::<Vec<PyParticle>>(),
        )
    }

    pub fn steps_taken_array<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.0.steps_taken_array().into_pyarray(py)
    }

    pub fn steps_stored_array<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.0.steps_stored_array().into_pyarray(py)
    }

    /// Utility method to avoid iterating on the Python side, which requires cloning.
    pub fn _initial_pzetas<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let vector = self
            .0
            .particles()
            .iter()
            .map(|p| p.initial_conditions().pzeta0().unwrap_or(f64::NAN))
            .collect::<Vec<f64>>();
        PyArray1::from_vec(py, vector)
    }

    /// Utility method to avoid iterating on the Python side, which requires cloning.
    pub fn _initial_energies<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let vector = self
            .0
            .particles()
            .iter()
            .map(|p| p.initial_energy().unwrap_or(f64::NAN))
            .collect::<Vec<f64>>();
        PyArray1::from_vec(py, vector)
    }

    /// Utility method to avoid iterating on the Python side, which requires cloning.
    pub fn _orbit_types<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        // maybe there is a better way to do this
        let vector = self
            .0
            .particles()
            .iter()
            .map(|p| format!("{:?}", p.orbit_type()))
            .collect::<Vec<String>>();
        PyList::new(py, vector.iter())
    }

    pub fn get_array<'py>(&self, py: Python<'py>, name: &str) -> Result<Bound<'py, PyArray1<f64>>> {
        let array = match name {
            "energy_array" => self.0.energy_array(),
            "energy_rsd_array" => self.0.energy_rsd_array(),
            "omega_theta_array" => self.0.omega_theta_array(),
            "omega_zeta_array" => self.0.omega_zeta_array(),
            "qkinetic_array" => self.0.qkinetic_array(),
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "Queue".into(),
                    attr: name.into(),
                });
            }
        };
        Ok(array.into_pyarray(py))
    }

    pub fn __len__(&mut self) -> usize {
        self.0.particle_count()
    }

    pub fn __getitem__(&self, index: usize) -> Option<PyParticle> {
        self.0
            .particles()
            .get(index)
            .map(|particle| PyParticle(particle.clone()))
    }
}

// ===============================================================================================

wrapper_debug_export!(PyQueue);
impl_py_repr!(PyQueue, pretty);
