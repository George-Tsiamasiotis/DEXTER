//! `dexter-equilibrium` newtypes, constructors and method exports.

use dexter::dexter_equilibrium::FluxCommute;
use dexter::dexter_equilibrium::{Current, LarCurrent, NcCurrent, NcCurrentBuilder};
use dexter::dexter_equilibrium::{Geometry, NcGeometry, NcGeometryBuilder};
use dexter::dexter_equilibrium::{
    NcQfactor, NcQfactorBuilder, ParabolicQfactor, Qfactor, UnityQfactor,
};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::PyEqError;
use crate::{
    py_debug_impl, py_eval1D, py_eval2D, py_export_getter, py_get_enum_string,
    py_get_netcdf_version, py_get_numpy1D, py_get_numpy2D, py_get_path, py_repr_impl,
};

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "NcGeometry", frozen)]
pub struct PyNcGeometry(NcGeometry);

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
py_export_getter!(PyNcGeometry, rwall, Option<f64>);
py_export_getter!(PyNcGeometry, psi_wall, Option<f64>);
py_export_getter!(PyNcGeometry, psip_wall, Option<f64>);
py_get_enum_string!(PyNcGeometry, psi_state);
py_get_enum_string!(PyNcGeometry, psip_state);
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
// ===============================================================================================

#[pyclass(name = "UnityQfactor", frozen)]
pub struct PyUnityQfactor(UnityQfactor);

#[pymethods]
impl PyUnityQfactor {
    #[new]
    #[pyo3(signature = ())]
    pub fn new() -> Self {
        Self(UnityQfactor::new())
    }
}

py_debug_impl!(PyUnityQfactor);
py_repr_impl!(PyUnityQfactor);
py_get_enum_string!(PyUnityQfactor, equilibrium_type);
py_eval1D!(PyUnityQfactor, q_of_psi);
py_eval1D!(PyUnityQfactor, q_of_psip);
py_eval1D!(PyUnityQfactor, psip_of_psi);
py_eval1D!(PyUnityQfactor, psi_of_psip);
py_eval1D!(PyUnityQfactor, dpsi_dpsip);
py_eval1D!(PyUnityQfactor, dpsip_dpsi);
py_eval1D!(PyUnityQfactor, iota_of_psi);
py_eval1D!(PyUnityQfactor, iota_of_psip);

// ===============================================================================================

#[pyclass(name = "ParabolicQfactor", frozen)]
pub struct PyParabolicQfactor(ParabolicQfactor);

#[pymethods]
impl PyParabolicQfactor {
    #[new]
    #[pyo3(signature = (qaxis, qwall, psi_wall))]
    pub fn new(qaxis: f64, qwall: f64, psi_wall: f64) -> Self {
        Self(ParabolicQfactor::new(qaxis, qwall, psi_wall))
    }
}

py_debug_impl!(PyParabolicQfactor);
py_repr_impl!(PyParabolicQfactor);
py_get_enum_string!(PyParabolicQfactor, equilibrium_type);
py_export_getter!(PyParabolicQfactor, qaxis, f64);
py_export_getter!(PyParabolicQfactor, qwall, f64);
py_export_getter!(PyParabolicQfactor, psi_wall, f64);
py_export_getter!(PyParabolicQfactor, psip_wall, f64);
py_eval1D!(PyParabolicQfactor, q_of_psi);
py_eval1D!(PyParabolicQfactor, q_of_psip);
py_eval1D!(PyParabolicQfactor, psip_of_psi);
py_eval1D!(PyParabolicQfactor, psi_of_psip);
py_eval1D!(PyParabolicQfactor, dpsi_dpsip);
py_eval1D!(PyParabolicQfactor, dpsip_dpsi);
py_eval1D!(PyParabolicQfactor, iota_of_psi);
py_eval1D!(PyParabolicQfactor, iota_of_psip);

// ===============================================================================================

#[pyclass(name = "NcQfactor", frozen)]
pub struct PyNcQfactor(NcQfactor);

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
py_export_getter!(PyNcQfactor, interp_type, String);
py_export_getter!(PyNcQfactor, qaxis, f64);
py_export_getter!(PyNcQfactor, qwall, f64);
py_export_getter!(PyNcQfactor, psi_wall, Option<f64>);
py_export_getter!(PyNcQfactor, psip_wall, Option<f64>);
py_get_enum_string!(PyNcQfactor, psi_state);
py_get_enum_string!(PyNcQfactor, psip_state);
py_get_numpy1D!(PyNcQfactor, psi_array);
py_get_numpy1D!(PyNcQfactor, psip_array);
py_get_numpy1D!(PyNcQfactor, q_array);
py_eval1D!(PyNcQfactor, q_of_psi);
py_eval1D!(PyNcQfactor, q_of_psip);
py_eval1D!(PyNcQfactor, psip_of_psi);
py_eval1D!(PyNcQfactor, psi_of_psip);
py_eval1D!(PyNcQfactor, dpsi_dpsip);
py_eval1D!(PyNcQfactor, dpsip_dpsi);
py_eval1D!(PyNcQfactor, iota_of_psi);
py_eval1D!(PyNcQfactor, iota_of_psip);

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "NcCurrent", frozen)]
pub struct PyNcCurrent(NcCurrent);

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
py_export_getter!(PyNcCurrent, psi_wall, Option<f64>);
py_export_getter!(PyNcCurrent, psip_wall, Option<f64>);
py_get_enum_string!(PyNcCurrent, psi_state);
py_get_enum_string!(PyNcCurrent, psip_state);
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

// ===============================================================================================

#[pyclass(name = "LarCurrent", frozen)]
pub struct PyLarCurrent(LarCurrent);

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
