//! Heap objects' Python wrappers.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use heap::{Heap, HeapInitialConditions};
use utils::{py_debug_impl, py_export_getter, py_get_numpy1D, py_repr_impl};

use crate::pyerrors::PyHeapError;
use crate::pylibrium::{PyBfield, PyCurrents, PyPerturbation, PyQfactor};
use crate::pyparticle::{PyMappingParameters, PyParticle};

#[pyclass(frozen, name = "HeapInitialConditions")]
pub struct PyHeapInitialConditions(pub HeapInitialConditions);

#[pymethods]
impl PyHeapInitialConditions {
    #[new]
    // Use Vec<f64> and let pyo3 do the rest.
    pub fn new_py(
        thetas: Vec<f64>,
        psips: Vec<f64>,
        rhos: Vec<f64>,
        zetas: Vec<f64>,
        mus: Vec<f64>,
    ) -> Result<Self, PyHeapError> {
        Ok(Self(HeapInitialConditions::build(
            &thetas, &psips, &rhos, &zetas, &mus,
        )?))
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}

py_debug_impl!(PyHeapInitialConditions);
py_repr_impl!(PyHeapInitialConditions);
py_get_numpy1D!(PyHeapInitialConditions, thetas);
py_get_numpy1D!(PyHeapInitialConditions, psips);
py_get_numpy1D!(PyHeapInitialConditions, rhos);
py_get_numpy1D!(PyHeapInitialConditions, zetas);
py_get_numpy1D!(PyHeapInitialConditions, mus);
py_export_getter!(PyHeapInitialConditions, len, usize);
py_export_getter!(PyHeapInitialConditions, is_empty, bool);

// ===============================================================================================

#[pyclass(name = "Heap")]
pub struct PyHeap(pub Heap);

#[pymethods]
impl PyHeap {
    #[new]
    pub fn new_py(initials: &PyHeapInitialConditions) -> Self {
        Self(Heap::new(&initials.0))
    }

    pub fn poincare(
        &mut self,
        qfactor: &PyQfactor,
        currents: &PyCurrents,
        bfield: &PyBfield,
        perturbation: &PyPerturbation,
        params: &PyMappingParameters,
    ) -> Result<(), PyHeapError> {
        Ok(self.0.poincare(
            &qfactor.0,
            &currents.0,
            &bfield.0,
            &perturbation.0,
            &params.0,
        )?)
    }

    /// Makes Heap indexable
    pub fn __getitem__(&self, index: usize) -> PyParticle {
        PyParticle(self.0.initials.particle_from_index(index))
    }

    pub fn __len__(&self) -> usize {
        self.0.len()
    }
}
