//! `dexter-equilibrium` newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{Current, NcCurrent, NcCurrentBuilder};
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use rsl_interpolation::Accelerator;

use crate::py_eval1D;
use crate::py_get_numpy1D;
use crate::pyerror::PyEqError;
use crate::{py_debug_impl, py_export_getter, py_get_path, py_len_impl, py_repr_impl};

// ===============================================================================================
// ===============================================================================================

#[pyclass(frozen, name = "NcCurrent")]
pub struct PyNcCurrent(pub NcCurrent);

#[pymethods]
impl PyNcCurrent {
    #[new]
    #[pyo3(signature = (path, typ))]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcCurrentBuilder::new(&path, typ);
        Ok(Self(builder.build()?))
    }
}

py_debug_impl!(PyNcCurrent);
py_repr_impl!(PyNcCurrent);
py_get_path!(PyNcCurrent);
py_len_impl!(PyNcCurrent);
py_export_getter!(PyNcCurrent, typ, String);
py_export_getter!(PyNcCurrent, psi_wall, Option<f64>);
py_export_getter!(PyNcCurrent, psip_wall, Option<f64>);
py_get_numpy1D!(PyNcCurrent, psi_array);
py_get_numpy1D!(PyNcCurrent, psip_array);
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
