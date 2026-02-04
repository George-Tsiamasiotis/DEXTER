//! `dexter-equilibrium` newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{Current, NcCurrent, NcCurrentBuilder};
use dexter::dexter_equilibrium::{Geometry, NcGeometry, NcGeometryBuilder};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::PyEqError;
use crate::{
    py_debug_impl, py_eval1D, py_eval2D, py_export_getter, py_get_enum_string, py_get_numpy1D,
    py_get_numpy2D, py_get_path, py_len_impl, py_repr_impl,
};

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "NcGeometry", frozen)]
pub struct PyNcGeometry(pub NcGeometry);

#[pymethods]
impl PyNcGeometry {
    #[new]
    #[pyo3(signature = (path, typ1d, typ2d))]
    pub fn new(path: &str, typ1d: &str, typ2d: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
        Ok(Self(builder.build()?))
    }

    #[getter]
    pub fn get_shape(&self) -> (usize, usize) {
        self.0.shape()
    }
}

py_debug_impl!(PyNcGeometry);
py_repr_impl!(PyNcGeometry);
py_get_path!(PyNcGeometry);
py_export_getter!(PyNcGeometry, typ1d, String);
py_export_getter!(PyNcGeometry, typ2d, String);
py_export_getter!(PyNcGeometry, baxis, f64);
py_export_getter!(PyNcGeometry, raxis, f64);
py_export_getter!(PyNcGeometry, zaxis, f64);
py_export_getter!(PyNcGeometry, rgeo, f64);
py_export_getter!(PyNcGeometry, rwall, Option<f64>);
py_get_enum_string!(PyNcGeometry, psi_state);
py_get_enum_string!(PyNcGeometry, psip_state);
py_export_getter!(PyNcGeometry, psi_wall, Option<f64>);
py_export_getter!(PyNcGeometry, psip_wall, Option<f64>);
py_get_numpy1D!(PyNcGeometry, psi_array);
py_get_numpy1D!(PyNcGeometry, psip_array);
py_get_numpy1D!(PyNcGeometry, theta_array);
py_get_numpy1D!(PyNcGeometry, r_array);
py_get_numpy2D!(PyNcGeometry, rlab_array);
py_get_numpy2D!(PyNcGeometry, zlab_array);
py_get_numpy2D!(PyNcGeometry, jacobian_array);
py_eval1D!(PyNcGeometry, psip_of_psi);
py_eval1D!(PyNcGeometry, psi_of_psip);
py_eval1D!(PyNcGeometry, r_of_psi);
py_eval1D!(PyNcGeometry, r_of_psip);
py_eval1D!(PyNcGeometry, psi_of_r);
py_eval1D!(PyNcGeometry, psip_of_r);
py_eval2D!(PyNcGeometry, rlab_of_psi);
py_eval2D!(PyNcGeometry, rlab_of_psip);
py_eval2D!(PyNcGeometry, zlab_of_psi);
py_eval2D!(PyNcGeometry, zlab_of_psip);
py_eval2D!(PyNcGeometry, jacobian_of_psi);
py_eval2D!(PyNcGeometry, jacobian_of_psip);

// ===============================================================================================

#[pyclass(name = "NcCurrent", frozen)]
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
py_get_enum_string!(PyNcCurrent, psi_state);
py_get_enum_string!(PyNcCurrent, psip_state);
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
