//! `dexter-equilibrium` newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{Bfield, LarBfield, NcBfield, NcBfieldBuilder};
use dexter::dexter_equilibrium::{Current, LarCurrent, NcCurrent, NcCurrentBuilder};
use dexter::dexter_equilibrium::{FluxCommute, FluxWall};
use dexter::dexter_equilibrium::{Geometry, LarGeometry, NcGeometry, NcGeometryBuilder};
use dexter::dexter_equilibrium::{
    NcQfactor, NcQfactorBuilder, ParabolicQfactor, Qfactor, UnityQfactor,
};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::PyEqError;
use crate::{
    py_debug_impl, py_eval1D, py_eval2D, py_export_getter, py_get_enum_string,
    py_get_netcdf_version, py_get_numpy1D, py_get_numpy1D_fallible, py_get_numpy2D, py_get_path,
    py_repr_impl,
};

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "_PyLarGeometry", subclass, frozen)]
pub struct PyLarGeometry(LarGeometry);

#[pymethods]
impl PyLarGeometry {
    #[new]
    #[pyo3(signature = (baxis, raxis, rwall))]
    pub fn new(baxis: f64, raxis: f64, rwall: f64) -> Self {
        Self(LarGeometry::new(baxis, raxis, rwall))
    }
}

py_debug_impl!(PyLarGeometry);
py_repr_impl!(PyLarGeometry);
py_get_enum_string!(PyLarGeometry, equilibrium_type);
py_export_getter!(PyLarGeometry, baxis, f64);
py_export_getter!(PyLarGeometry, raxis, f64);
py_export_getter!(PyLarGeometry, rwall, f64);
py_export_getter!(PyLarGeometry, psi_wall, f64);
py_get_numpy1D!(PyLarGeometry, rlab_wall);
py_get_numpy1D!(PyLarGeometry, zlab_wall);
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

// ===============================================================================================

#[pyclass(name = "_PyNcGeometry", subclass, frozen)]
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
py_export_getter!(PyNcGeometry, rwall, f64);
py_export_getter!(PyNcGeometry, psi_wall, Option<f64>);
py_export_getter!(PyNcGeometry, psip_wall, Option<f64>);
py_get_enum_string!(PyNcGeometry, psi_state);
py_get_enum_string!(PyNcGeometry, psip_state);
py_get_numpy1D_fallible!(PyNcGeometry, psi_array);
py_get_numpy1D_fallible!(PyNcGeometry, psip_array);
py_get_numpy1D!(PyNcGeometry, theta_array);
py_get_numpy1D!(PyNcGeometry, r_array);
py_get_numpy2D!(PyNcGeometry, rlab_array);
py_get_numpy2D!(PyNcGeometry, zlab_array);
py_get_numpy2D!(PyNcGeometry, jacobian_array);
py_get_numpy1D!(PyNcGeometry, rlab_wall);
py_get_numpy1D!(PyNcGeometry, zlab_wall);
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

#[pyclass(name = "_PyUnityQfactor", subclass, frozen)]
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

#[pyclass(name = "_PyParabolicQfactor", subclass, frozen)]
pub struct PyParabolicQfactor(ParabolicQfactor);

#[pymethods]
impl PyParabolicQfactor {
    #[new]
    #[pyo3(signature = (qaxis, qwall, flux_wall))]
    pub fn new<'py>(qaxis: f64, qwall: f64, flux_wall: Bound<'py, PyTuple>) -> Self {
        let e = "Input type must be a tuple of the form (str, float)";
        let flux: String = flux_wall.get_item(0).expect(e).extract().expect(e);
        let value: f64 = flux_wall.get_item(1).expect(e).extract().expect(e);
        let flux_wall: FluxWall;
        match flux.to_lowercase().as_str() {
            "toroidal" => flux_wall = FluxWall::Toroidal(value),
            "poloidal" => flux_wall = FluxWall::Poloidal(value),
            _ => panic!("Flux coordinate must be either 'Toroidal' or 'Poloidal'"),
        }

        Self(ParabolicQfactor::new(qaxis, qwall, flux_wall))
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

#[pyclass(name = "_PyNcQfactor", subclass, frozen)]
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
py_get_numpy1D_fallible!(PyNcQfactor, psi_array);
py_get_numpy1D_fallible!(PyNcQfactor, psip_array);
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

#[pyclass(name = "_PyLarCurrent", subclass, frozen)]
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

// ===============================================================================================

#[pyclass(name = "_PyNcCurrent", subclass, frozen)]
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

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "_PyLarBfield", subclass, frozen)]
pub struct PyLarBfield(LarBfield);

#[pymethods]
impl PyLarBfield {
    #[new]
    #[pyo3(signature = ())]
    pub fn new() -> Self {
        Self(LarBfield::new())
    }
}

py_debug_impl!(PyLarBfield);
py_repr_impl!(PyLarBfield);
py_get_enum_string!(PyLarBfield, equilibrium_type);
py_eval2D!(PyLarBfield, b_of_psi);
py_eval2D!(PyLarBfield, b_of_psip);
py_eval2D!(PyLarBfield, db_dpsi);
py_eval2D!(PyLarBfield, db_dpsip);
py_eval2D!(PyLarBfield, db_of_psi_dtheta);
py_eval2D!(PyLarBfield, db_of_psip_dtheta);

// ===============================================================================================

#[pyclass(name = "_PyNcBfield", subclass, frozen)]
pub struct PyNcBfield(NcBfield);

#[pymethods]
impl PyNcBfield {
    #[new]
    #[pyo3(signature = (path, interp_type))]
    pub fn new(path: &str, interp_type: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcBfieldBuilder::new(&path, interp_type);
        Ok(Self(builder.build()?))
    }

    #[getter]
    pub fn get_shape(&self) -> (usize, usize) {
        self.0.shape()
    }
}

py_debug_impl!(PyNcBfield);
py_repr_impl!(PyNcBfield);
py_get_path!(PyNcBfield);
py_get_netcdf_version!(PyNcBfield);
py_get_enum_string!(PyNcBfield, equilibrium_type);
py_export_getter!(PyNcBfield, interp_type, String);
py_export_getter!(PyNcBfield, baxis, f64);
py_export_getter!(PyNcBfield, psi_wall, Option<f64>);
py_export_getter!(PyNcBfield, psip_wall, Option<f64>);
py_get_enum_string!(PyNcBfield, psi_state);
py_get_enum_string!(PyNcBfield, psip_state);
py_get_numpy1D_fallible!(PyNcBfield, psi_array);
py_get_numpy1D_fallible!(PyNcBfield, psip_array);
py_get_numpy1D!(PyNcBfield, theta_array);
py_get_numpy2D!(PyNcBfield, b_array);
py_eval2D!(PyNcBfield, b_of_psi);
py_eval2D!(PyNcBfield, b_of_psip);
py_eval2D!(PyNcBfield, db_dpsi);
py_eval2D!(PyNcBfield, db_dpsip);
py_eval2D!(PyNcBfield, db_of_psi_dtheta);
py_eval2D!(PyNcBfield, db_of_psip_dtheta);
