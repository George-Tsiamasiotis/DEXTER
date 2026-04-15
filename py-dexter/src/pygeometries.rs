//! `dexter-equilibrium` geometries' newtypes, constructors and method exports.

use dexter::dexter_equilibrium::FluxCommute;
use dexter::dexter_equilibrium::{Geometry, LarGeometry, NcGeometry, NcGeometryBuilder};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::{PyEqError, PyEvalError};
use crate::{
    py_debug_impl, py_eval1D, py_eval2D, py_export_getter, py_get_enum_string,
    py_get_netcdf_version, py_get_numpy1D, py_get_numpy1D_fallible, py_get_numpy2D, py_get_path,
    py_repr_impl,
};

// ===============================================================================================

#[pyclass(name = "_PyLarGeometry", frozen)]
pub struct PyLarGeometry(pub(crate) LarGeometry);

#[pymethods]
impl PyLarGeometry {
    #[new]
    #[pyo3(signature = (baxis, raxis, rlast))]
    pub fn new(baxis: f64, raxis: f64, rlast: f64) -> Self {
        Self(LarGeometry::new(baxis, raxis, rlast))
    }
}

py_debug_impl!(PyLarGeometry);
py_repr_impl!(PyLarGeometry);
py_get_enum_string!(PyLarGeometry, equilibrium_type);
py_export_getter!(PyLarGeometry, baxis, f64);
py_export_getter!(PyLarGeometry, raxis, f64);
py_export_getter!(PyLarGeometry, zaxis, f64);
py_export_getter!(PyLarGeometry, rgeo, f64);
py_export_getter!(PyLarGeometry, rlast, f64);
py_export_getter!(PyLarGeometry, psi_last, f64);
py_eval1D!(PyLarGeometry, r_of_psi);
py_eval1D!(PyLarGeometry, r_of_psip);
py_eval1D!(PyLarGeometry, psi_of_r);
py_eval1D!(PyLarGeometry, psip_of_r);
py_eval2D!(PyLarGeometry, rlab_of_psi);
py_eval2D!(PyLarGeometry, rlab_of_psip);
py_eval2D!(PyLarGeometry, zlab_of_psi);
py_eval2D!(PyLarGeometry, zlab_of_psip);
py_eval2D!(PyLarGeometry, jacobian_of_psi);
py_eval2D!(PyLarGeometry, jacobian_of_psip);
py_get_numpy1D!(PyLarGeometry, rlab_last);
py_get_numpy1D!(PyLarGeometry, zlab_last);

// ===============================================================================================

#[pyclass(name = "_PyNcGeometry", frozen)]
pub struct PyNcGeometry(pub(crate) NcGeometry);

#[pymethods]
impl PyNcGeometry {
    #[new]
    #[pyo3(signature = (path, interp1d_type, interp2d_type))]
    pub fn new(path: &str, interp1d_type: &str, interp2d_type: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type);
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
py_get_netcdf_version!(PyNcGeometry);
py_get_enum_string!(PyNcGeometry, equilibrium_type);
py_export_getter!(PyNcGeometry, interp1d_type, String);
py_export_getter!(PyNcGeometry, interp2d_type, String);
py_export_getter!(PyNcGeometry, baxis, f64);
py_export_getter!(PyNcGeometry, raxis, f64);
py_export_getter!(PyNcGeometry, zaxis, f64);
py_export_getter!(PyNcGeometry, rgeo, f64);
py_export_getter!(PyNcGeometry, rlast, f64);
py_export_getter!(PyNcGeometry, psi_last, Option<f64>);
py_export_getter!(PyNcGeometry, psip_last, Option<f64>);
py_get_enum_string!(PyNcGeometry, psi_state);
py_get_enum_string!(PyNcGeometry, psip_state);
py_get_numpy1D_fallible!(PyNcGeometry, psi_array);
py_get_numpy1D_fallible!(PyNcGeometry, psip_array);
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
py_get_numpy1D!(PyNcGeometry, rlab_last);
py_get_numpy1D!(PyNcGeometry, zlab_last);
