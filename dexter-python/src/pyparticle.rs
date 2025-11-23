//! Particle objects' Python wrappers.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use particle::{Evolution, Frequencies, InitialConditions};
use utils::{
    py_debug_impl, py_export_getter, py_get_numpy1D, py_get_primitive_field, py_repr_impl,
};

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
