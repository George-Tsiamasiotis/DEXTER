//! Particle objects' Python wrappers.

use pyo3::prelude::*;

use particle::InitialConditions;
use utils::{py_debug_impl, py_get_primitive_field, py_repr_impl};

#[pyclass(frozen, name = "InitialConditions")]
pub struct PyInitialConditions(InitialConditions);

#[pymethods]
impl PyInitialConditions {
    /// Creates a new PyQFactor wrapper object.
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
