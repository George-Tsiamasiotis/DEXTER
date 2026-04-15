//! `dexter-equilibrium` bfields' newtypes, constructors and method exports.

use dexter::dexter_equilibrium::{Bfield, LarBfield, NcBfield, NcBfieldBuilder};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyInt;
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::{PyEqError, PyEvalError};
use crate::{
    py_debug_impl, py_eval2D, py_export_getter, py_get_enum_string, py_get_netcdf_version,
    py_get_numpy1D, py_get_numpy1D_fallible, py_get_numpy2D, py_get_path, py_repr_impl,
};

// ===============================================================================================

#[pyclass(name = "_PyLarBfield", frozen)]
pub struct PyLarBfield(pub(crate) LarBfield);

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

#[pyclass(name = "_PyNcBfield", frozen)]
pub struct PyNcBfield(pub(crate) NcBfield);

#[pymethods]
impl PyNcBfield {
    #[new]
    #[pyo3(signature = (path, interp_type, *, padding))]
    pub fn new<'py>(
        path: &str,
        interp_type: &str,
        padding: Option<Bound<'py, PyInt>>,
    ) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let _padding = match padding {
            Some(py_int) => py_int.extract::<usize>()?,
            None => 15,
        };
        let builder = NcBfieldBuilder::new(&path, interp_type).with_padding(_padding);
        Ok(Self(builder.build()?))
    }

    #[getter]
    pub fn get_shape(&self) -> (usize, usize) {
        self.0.shape()
    }

    #[getter]
    pub fn get_shape_padded(&self) -> (usize, usize) {
        self.0.shape_padded()
    }
}

py_debug_impl!(PyNcBfield);
py_repr_impl!(PyNcBfield);
py_get_path!(PyNcBfield);
py_get_netcdf_version!(PyNcBfield);
py_get_enum_string!(PyNcBfield, equilibrium_type);
py_export_getter!(PyNcBfield, interp_type, String);
py_export_getter!(PyNcBfield, baxis, f64);
py_export_getter!(PyNcBfield, padding, usize);
py_export_getter!(PyNcBfield, psi_last, Option<f64>);
py_export_getter!(PyNcBfield, psip_last, Option<f64>);
py_get_enum_string!(PyNcBfield, psi_state);
py_get_enum_string!(PyNcBfield, psip_state);
py_get_numpy1D_fallible!(PyNcBfield, psi_array);
py_get_numpy1D_fallible!(PyNcBfield, psip_array);
py_get_numpy1D!(PyNcBfield, theta_array);
py_get_numpy2D!(PyNcBfield, b_array);
py_get_numpy1D!(PyNcBfield, theta_array_padded);
py_get_numpy2D!(PyNcBfield, b_array_padded);
py_eval2D!(PyNcBfield, b_of_psi);
py_eval2D!(PyNcBfield, b_of_psip);
py_eval2D!(PyNcBfield, db_dpsi);
py_eval2D!(PyNcBfield, db_dpsip);
py_eval2D!(PyNcBfield, db_of_psi_dtheta);
py_eval2D!(PyNcBfield, db_of_psip_dtheta);
