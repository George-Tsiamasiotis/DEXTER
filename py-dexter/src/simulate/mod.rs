mod energy;
mod initial;
mod particle;
mod queue;

pub use energy::*;
pub use initial::*;
pub use particle::*;
pub use queue::*;

use dexter::dexter_simulate::{IntersectParams, Intersection, SolverParams, SteppingMethod};
use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::*;

// ===============================================================================================

#[pyclass(name = "_PySteppingMethod", frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PySteppingMethod(pub SteppingMethod);

#[pymethods]
impl PySteppingMethod {
    #[classmethod]
    pub fn energy_adaptive_step(_: &Bound<'_, PyType>, rel_tol: f64, abs_tol: f64) -> Self {
        Self(SteppingMethod::EnergyAdaptiveStep { rel_tol, abs_tol })
    }

    #[classmethod]
    pub fn error_adaptive_step(_: &Bound<'_, PyType>, rel_tol: f64, abs_tol: f64) -> Self {
        Self(SteppingMethod::ErrorAdaptiveStep { rel_tol, abs_tol })
    }

    #[classmethod]
    pub fn fixed_step(_: &Bound<'_, PyType>, step: f64) -> Self {
        Self(SteppingMethod::FixedStep(step))
    }
}

/// This type is only to be used internally when calling integration routines.
#[pyclass(name = "_PySolverParams", frozen, immutable_type)]
pub struct PySolverParams(SolverParams);

#[pymethods]
impl PySolverParams {
    #[new]
    pub fn new<'py>(
        method: Option<PySteppingMethod>,
        max_steps: Option<usize>,
        first_step: Option<f64>,
        safety_factor: Option<f64>,
    ) -> PyResult<Self> {
        let mut solver_params = SolverParams::default();
        method.inspect(|v| solver_params.method = v.clone().0);
        max_steps.inspect(|v| solver_params.max_steps = *v);
        first_step.inspect(|v| solver_params.first_step = *v);
        safety_factor.inspect(|v| solver_params.safety_factor = *v);

        Ok(Self(solver_params))
    }
}

// ===============================================================================================

#[pyclass(name = "_PyIntersectParams", frozen, immutable_type)]
pub struct PyIntersectParams(pub(crate) IntersectParams);

#[pymethods]
impl PyIntersectParams {
    #[new]
    pub fn new<'py>(intersection: String, angle: f64, turns: usize) -> PyResult<Self> {
        let intersection = match intersection.to_lowercase().as_str() {
            "consttheta" => Intersection::ConstTheta,
            "constzeta" => Intersection::ConstZeta,
            _ => return Err(PyErr::from(DexterError::InvalidIntersection)),
        };
        let intersect_params = IntersectParams::new(intersection, angle, turns);

        Ok(Self(intersect_params))
    }
}

// ===============================================================================================

wrapper_debug_export!(PySteppingMethod);
wrapper_debug_export!(PySolverParams);
wrapper_debug_export!(PyIntersectParams);
impl_py_repr!(PySteppingMethod, pretty);
impl_py_repr!(PySolverParams, pretty);
impl_py_repr!(PyIntersectParams, pretty);
