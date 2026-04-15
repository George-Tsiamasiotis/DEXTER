//! `dexter-equilibrium` harmonics' newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{
    CosHarmonic, Harmonic, NcHarmonic, NcHarmonicBuilder, PhaseMethod,
};
use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::pyerror::{PyEqError, PyEvalError};
use crate::pylibrium_misc::PyLastClosedFluxSurface;
use crate::{
    py_debug_impl, py_eval_harmonic, py_export_getter, py_get_enum_string, py_get_netcdf_version,
    py_get_numpy1D, py_get_numpy1D_fallible, py_get_path, py_repr_impl,
};

// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyCosHarmonic", frozen)]
pub struct PyCosHarmonic(pub(crate) CosHarmonic);

#[pymethods]
impl PyCosHarmonic {
    #[new]
    #[pyo3(signature = (epsilon, lcfs, m, n, phase))]
    pub fn new(epsilon: f64, lcfs: &PyLastClosedFluxSurface, m: i64, n: i64, phase: f64) -> Self {
        Self(CosHarmonic::new(epsilon, lcfs.0, m, n, phase))
    }

    #[getter]
    pub fn get_lcfs(&self) -> PyLastClosedFluxSurface {
        PyLastClosedFluxSurface(self.0.lcfs())
    }
}

py_debug_impl!(PyCosHarmonic);
py_repr_impl!(PyCosHarmonic);
py_get_enum_string!(PyCosHarmonic, equilibrium_type);
py_export_getter!(PyCosHarmonic, epsilon, f64);
py_export_getter!(PyCosHarmonic, psi_last, Option<f64>);
py_export_getter!(PyCosHarmonic, psip_last, Option<f64>);
py_export_getter!(PyCosHarmonic, phase, f64);
py_export_getter!(PyCosHarmonic, m, i64);
py_export_getter!(PyCosHarmonic, n, i64);
py_eval_harmonic!(PyCosHarmonic, alpha_of_psi);
py_eval_harmonic!(PyCosHarmonic, alpha_of_psip);
py_eval_harmonic!(PyCosHarmonic, phase_of_psi);
py_eval_harmonic!(PyCosHarmonic, phase_of_psip);
py_eval_harmonic!(PyCosHarmonic, h_of_psi);
py_eval_harmonic!(PyCosHarmonic, h_of_psip);
py_eval_harmonic!(PyCosHarmonic, dh_dpsi);
py_eval_harmonic!(PyCosHarmonic, dh_dpsip);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dtheta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dtheta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dzeta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dzeta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dt);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dt);

// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyNcHarmonic", frozen)]
pub struct PyNcHarmonic(pub(crate) NcHarmonic);

#[pymethods]
impl PyNcHarmonic {
    #[new]
    #[pyo3(signature = (path, interp_type, m, n, phase_method))]
    pub fn new<'py>(
        path: &str,
        interp_type: &str,
        m: i64,
        n: i64,
        phase_method: Option<Bound<'py, PyAny>>,
    ) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let method: PhaseMethod = Self::resolve_phase_method(phase_method)?;
        let builder = NcHarmonicBuilder::new(&path, interp_type, m, n).with_phase_method(method);
        Ok(Self(builder.build()?))
    }
}

impl PyNcHarmonic {
    /// Attempts to create a valid [`PhaseMethod`] object from the corresponding [`PyAny`]
    /// Python argument.
    ///
    /// Returns early if a valid string is found that matches one of the simple variants,
    /// otherwise tries to unpack the [`PyAny`] object into a ("Custom", f64) tuple and cast it
    /// into a PhaseMethod::Custom.
    ///
    /// Returns an Error if both string matching and casting fail.
    fn resolve_phase_method<'py>(arg: Option<Bound<'py, PyAny>>) -> PyResult<PhaseMethod> {
        use PhaseMethod::*;

        if arg.is_none() {
            return Ok(PhaseMethod::Average);
        }
        let arg = arg.unwrap();

        match arg.to_string().to_lowercase().as_str() {
            "zero" => return Ok(Zero),
            "average" => return Ok(Average),
            "resonance" => return Ok(Resonance),
            "interpolation" => return Ok(Interpolation),
            _ => (),
        };

        let tuple = arg.cast::<PyTuple>()?;
        let string: String = tuple.get_item(0)?.extract::<String>()?.to_lowercase();
        let value: f64 = tuple.get_item(1)?.extract::<f64>()?;
        match string.as_str() {
            "custom" if value.is_finite() => Ok(PhaseMethod::Custom(value)),
            _ => Err(PyErr::new::<PyTypeError, _>("Invalid 'phase_method'")),
        }
    }
}

py_debug_impl!(PyNcHarmonic);
py_repr_impl!(PyNcHarmonic);
py_get_path!(PyNcHarmonic);
py_get_netcdf_version!(PyNcHarmonic);
py_get_enum_string!(PyNcHarmonic, equilibrium_type);
py_export_getter!(PyNcHarmonic, interp_type, String);
py_export_getter!(PyNcHarmonic, m, i64);
py_export_getter!(PyNcHarmonic, n, i64);
py_get_enum_string!(PyNcHarmonic, phase_method);
py_export_getter!(PyNcHarmonic, phase_average, Option<f64>);
py_export_getter!(PyNcHarmonic, psi_phase_resonance, Option<f64>);
py_export_getter!(PyNcHarmonic, psip_phase_resonance, Option<f64>);
py_export_getter!(PyNcHarmonic, psi_last, Option<f64>);
py_export_getter!(PyNcHarmonic, psip_last, Option<f64>);
py_get_enum_string!(PyNcHarmonic, psi_state);
py_get_enum_string!(PyNcHarmonic, psip_state);
py_get_numpy1D_fallible!(PyNcHarmonic, psi_array);
py_get_numpy1D_fallible!(PyNcHarmonic, psip_array);
py_get_numpy1D!(PyNcHarmonic, alpha_array);
py_get_numpy1D!(PyNcHarmonic, phase_array);
py_eval_harmonic!(PyNcHarmonic, alpha_of_psi);
py_eval_harmonic!(PyNcHarmonic, alpha_of_psip);
py_eval_harmonic!(PyNcHarmonic, phase_of_psi);
py_eval_harmonic!(PyNcHarmonic, phase_of_psip);
py_eval_harmonic!(PyNcHarmonic, h_of_psi);
py_eval_harmonic!(PyNcHarmonic, h_of_psip);
py_eval_harmonic!(PyNcHarmonic, dh_dpsi);
py_eval_harmonic!(PyNcHarmonic, dh_dpsip);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dtheta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dtheta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dzeta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dzeta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dt);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dt);
