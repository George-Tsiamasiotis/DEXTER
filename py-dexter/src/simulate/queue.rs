//! Defines the wrapper around [`dexter::dexter_simulate::Queue`].

use dexter::dexter_machine::*;
use dexter::dexter_simulate::*;

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
        solver_params: &PySolverParams,
    ) {
        let machine = MachineBuilder::new(qfactor.inner(), current.inner(), bfield.inner())
            .with_perturbation(perturbation.inner())
            .build();
        self.0.close(machine, periods, &solver_params.0);
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

// ===============================================================================================

wrapper_debug_export!(PyQueue);
impl_py_repr!(PyQueue, pretty);
