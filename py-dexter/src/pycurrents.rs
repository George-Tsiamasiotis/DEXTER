//! `dexter-equilibrium` currents' newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{Current, LarCurrent, NcCurrent, NcCurrentBuilder};
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use rsl_interpolation::Accelerator;

use crate::pyerror::{PyEqError, PyEvalError};
use crate::{
    py_debug_impl, py_eval1D, py_export_getter, py_get_enum_string, py_get_netcdf_version,
    py_get_numpy1D, py_get_numpy1D_fallible, py_get_path, py_repr_impl,
};

// ===============================================================================================

#[pyclass(name = "_PyLarCurrent", frozen)]
pub struct PyLarCurrent(pub(crate) LarCurrent);

#[pymethods]
impl PyLarCurrent {
    #[new]
    #[pyo3(signature = ())]
    pub fn new() -> Self {
        Self(LarCurrent::new())
    }
}

py_debug_impl!(PyLarCurrent);
py_repr_impl!(PyLarCurrent);
py_get_enum_string!(PyLarCurrent, equilibrium_type);
py_eval1D!(PyLarCurrent, g_of_psi);
py_eval1D!(PyLarCurrent, g_of_psip);
py_eval1D!(PyLarCurrent, dg_dpsi);
py_eval1D!(PyLarCurrent, dg_dpsip);
py_eval1D!(PyLarCurrent, i_of_psi);
py_eval1D!(PyLarCurrent, i_of_psip);
py_eval1D!(PyLarCurrent, di_dpsi);
py_eval1D!(PyLarCurrent, di_dpsip);

// ===============================================================================================

#[pyclass(name = "_PyNcCurrent", frozen)]
pub struct PyNcCurrent(pub(crate) NcCurrent);

#[pymethods]
impl PyNcCurrent {
    #[new]
    #[pyo3(signature = (path, interp_type))]
    pub fn new(path: &str, interp_type: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcCurrentBuilder::new(&path, interp_type);
        Ok(Self(builder.build()?))
    }
}

py_debug_impl!(PyNcCurrent);
py_repr_impl!(PyNcCurrent);
py_get_path!(PyNcCurrent);
py_get_netcdf_version!(PyNcCurrent);
py_get_enum_string!(PyNcCurrent, equilibrium_type);
py_export_getter!(PyNcCurrent, interp_type, String);
py_export_getter!(PyNcCurrent, psi_last, Option<f64>);
py_export_getter!(PyNcCurrent, psip_last, Option<f64>);
py_get_enum_string!(PyNcCurrent, psi_state);
py_get_enum_string!(PyNcCurrent, psip_state);
py_get_numpy1D_fallible!(PyNcCurrent, psi_array);
py_get_numpy1D_fallible!(PyNcCurrent, psip_array);
py_get_numpy1D!(PyNcCurrent, g_array);
py_get_numpy1D!(PyNcCurrent, i_array);
py_eval1D!(PyNcCurrent, g_of_psi);
py_eval1D!(PyNcCurrent, g_of_psip);
py_eval1D!(PyNcCurrent, dg_dpsi);
py_eval1D!(PyNcCurrent, dg_dpsip);
py_eval1D!(PyNcCurrent, i_of_psi);
py_eval1D!(PyNcCurrent, i_of_psip);
py_eval1D!(PyNcCurrent, di_dpsi);
py_eval1D!(PyNcCurrent, di_dpsip);
