//! `dexter-equilibrium` qfactors' newtypes, constructors and method exports.

use dexter::dexter_equilibrium::FluxCommute;
use dexter::dexter_equilibrium::{
    NcQfactor, NcQfactorBuilder, ParabolicQfactor, Qfactor, UnityQfactor,
};
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use rsl_interpolation::Accelerator;

use crate::pyerror::{PyEqError, PyEvalError};
use crate::pylibrium_misc::PyLastClosedFluxSurface;
use crate::{
    py_debug_impl, py_eval1D, py_export_getter, py_get_enum_string, py_get_netcdf_version,
    py_get_numpy1D, py_get_numpy1D_fallible, py_get_path, py_repr_impl,
};

// ===============================================================================================

#[pyclass(name = "_PyUnityQfactor", frozen)]
pub struct PyUnityQfactor(pub(crate) UnityQfactor);

#[pymethods]
impl PyUnityQfactor {
    #[new]
    #[pyo3(signature = (lcfs))]
    pub fn new(lcfs: &PyLastClosedFluxSurface) -> Self {
        Self(UnityQfactor::new(lcfs.0))
    }
}

py_debug_impl!(PyUnityQfactor);
py_repr_impl!(PyUnityQfactor);
py_get_enum_string!(PyUnityQfactor, equilibrium_type);
py_export_getter!(PyUnityQfactor, psi_last, f64);
py_export_getter!(PyUnityQfactor, psip_last, f64);
py_export_getter!(PyUnityQfactor, qlast, f64);
py_export_getter!(PyUnityQfactor, qaxis, f64);
py_eval1D!(PyUnityQfactor, psip_of_psi);
py_eval1D!(PyUnityQfactor, psi_of_psip);
py_eval1D!(PyUnityQfactor, q_of_psi);
py_eval1D!(PyUnityQfactor, q_of_psip);
py_eval1D!(PyUnityQfactor, dpsi_dpsip);
py_eval1D!(PyUnityQfactor, dpsip_dpsi);
py_eval1D!(PyUnityQfactor, iota_of_psi);
py_eval1D!(PyUnityQfactor, iota_of_psip);
py_eval1D!(PyUnityQfactor, psi_of_q);
py_eval1D!(PyUnityQfactor, psip_of_q);

// ===============================================================================================

#[pyclass(name = "_PyParabolicQfactor", frozen)]
pub struct PyParabolicQfactor(pub(crate) ParabolicQfactor);

#[pymethods]
impl PyParabolicQfactor {
    #[new]
    #[pyo3(signature = (qaxis, qlast, lcfs))]
    pub fn new<'py>(qaxis: f64, qlast: f64, lcfs: &PyLastClosedFluxSurface) -> PyResult<Self> {
        Ok(Self(ParabolicQfactor::new(qaxis, qlast, lcfs.0)))
    }
}

py_debug_impl!(PyParabolicQfactor);
py_repr_impl!(PyParabolicQfactor);
py_get_enum_string!(PyParabolicQfactor, equilibrium_type);
py_export_getter!(PyParabolicQfactor, psi_last, f64);
py_export_getter!(PyParabolicQfactor, psip_last, f64);
py_export_getter!(PyParabolicQfactor, qlast, f64);
py_export_getter!(PyParabolicQfactor, qaxis, f64);
py_eval1D!(PyParabolicQfactor, psip_of_psi);
py_eval1D!(PyParabolicQfactor, psi_of_psip);
py_eval1D!(PyParabolicQfactor, q_of_psi);
py_eval1D!(PyParabolicQfactor, q_of_psip);
py_eval1D!(PyParabolicQfactor, dpsi_dpsip);
py_eval1D!(PyParabolicQfactor, dpsip_dpsi);
py_eval1D!(PyParabolicQfactor, iota_of_psi);
py_eval1D!(PyParabolicQfactor, iota_of_psip);
py_eval1D!(PyParabolicQfactor, psi_of_q);
py_eval1D!(PyParabolicQfactor, psip_of_q);

// ===============================================================================================

#[pyclass(name = "_PyNcQfactor", frozen)]
pub struct PyNcQfactor(pub(crate) NcQfactor);

#[pymethods]
impl PyNcQfactor {
    #[new]
    #[pyo3(signature = (path, interp_type))]
    pub fn new(path: &str, interp_type: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcQfactorBuilder::new(&path, interp_type);
        Ok(Self(builder.build()?))
    }
}

py_debug_impl!(PyNcQfactor);
py_repr_impl!(PyNcQfactor);
py_get_path!(PyNcQfactor);
py_get_netcdf_version!(PyNcQfactor);
py_get_enum_string!(PyNcQfactor, equilibrium_type);
py_export_getter!(PyNcQfactor, psi_last, f64);
py_export_getter!(PyNcQfactor, psip_last, f64);
py_export_getter!(PyNcQfactor, qlast, f64);
py_export_getter!(PyNcQfactor, qaxis, f64);
py_export_getter!(PyNcQfactor, interp_type, String);
py_get_enum_string!(PyNcQfactor, psi_state);
py_get_enum_string!(PyNcQfactor, psip_state);
py_get_numpy1D_fallible!(PyNcQfactor, psi_array);
py_get_numpy1D_fallible!(PyNcQfactor, psip_array);
py_get_numpy1D!(PyNcQfactor, q_array);
py_eval1D!(PyNcQfactor, psip_of_psi);
py_eval1D!(PyNcQfactor, psi_of_psip);
py_eval1D!(PyNcQfactor, q_of_psi);
py_eval1D!(PyNcQfactor, q_of_psip);
py_eval1D!(PyNcQfactor, dpsi_dpsip);
py_eval1D!(PyNcQfactor, dpsip_dpsi);
py_eval1D!(PyNcQfactor, iota_of_psi);
py_eval1D!(PyNcQfactor, iota_of_psip);
py_eval1D!(PyNcQfactor, psi_of_q);
py_eval1D!(PyNcQfactor, psip_of_q);
