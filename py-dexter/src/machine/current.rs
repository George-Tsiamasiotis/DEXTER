//! Defines the PyCurrent enum that holds one of the Current objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyLarCurrent(Arc<LarCurrent>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcCurrent(Arc<NcCurrent>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyCurrent", frozen, immutable_type)]
pub enum PyCurrent {
    Lar(PyLarCurrent),
    Nc(PyNcCurrent),
}

#[pymethods] // Builders
impl PyCurrent {
    #[classmethod]
    pub fn build_lar<'py>(_: Bound<'py, PyType>) -> Result<Self> {
        Ok(Self::Lar(PyLarCurrent(Arc::new(LarCurrent::new()))))
    }

    #[classmethod]
    pub fn build_nc<'py>(_: Bound<'py, PyType>, path: String, interp_type: String) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_1d_type(interp_type)?;
        let builder = NcCurrentBuilder::new(&path, typ);
        let current = builder.build()?;
        Ok(Self::Nc(PyNcCurrent(Arc::new(current))))
    }
}

/// References to the trait object and variants
impl PyCurrent {
    pub fn inner(&self) -> &dyn Current {
        match self {
            PyCurrent::Lar(current) => current.0.as_ref(),
            PyCurrent::Nc(current) => current.0.as_ref(),
        }
    }

    pub fn lar(&self) -> Result<&LarCurrent> {
        match self {
            Self::Lar(current) => Ok(current.0.as_ref()),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Current".into(),
                inner: "LarCurrent".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcCurrent> {
        match self {
            Self::Nc(current) => Ok(&current.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Current".into(),
                inner: "NcCurrent".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // MachineObject Trait
impl PyCurrent {
    #[getter]
    pub fn machine_type(&self) -> String {
        format!("{:?}", self.inner().machine_type())
    }

    #[getter]
    pub fn psi_state(&self) -> String {
        format!("{:?}", self.inner().psi_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.inner().psip_state())
    }
}

#[pymethods] // Current Trait
impl PyCurrent {
    pub fn eval_g(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_g(flux, &mut Accelerator::new())?)
    }

    pub fn eval_i(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_i(flux, &mut Accelerator::new())?)
    }

    pub fn eval_g_deriv(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_g_deriv(flux, &mut Accelerator::new())?)
    }

    pub fn eval_i_deriv(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_i_deriv(flux, &mut Accelerator::new())?)
    }
}

// ===============================================================================================

// #[pymethods] // Lar
impl PyCurrent {}

#[pymethods] // Nc
impl PyCurrent {
    #[getter]
    pub fn path(&self) -> Result<String> {
        Ok(self.nc()?.path().to_str().unwrap_or_default().to_string())
    }

    #[getter]
    pub fn netcdf_version(&self) -> Result<String> {
        Ok(self.nc()?.netcdf_version().to_string())
    }

    #[getter]
    pub fn interp_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp_type()))
    }

    pub fn get_array<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        let current = self.nc()?;
        match name {
            "g_array" => Ok(Some(current.g_array().into_pyarray(py))),
            "i_array" => Ok(Some(current.i_array().into_pyarray(py))),
            "psi_array" => match current.psi_array() {
                Some(array) => Ok(Some(array.into_pyarray(py))),
                None => Ok(None),
            },
            "psip_array" => match current.psip_array() {
                Some(array) => Ok(Some(array.into_pyarray(py))),
                None => Ok(None),
            },
            _ => Err(DexterError::AttributeError {
                obj: "NcCurrent".into(),
                attr: name.into(),
            }),
        }
    }
}

// ===============================================================================================

wrapper_debug_export!(PyLarCurrent);
wrapper_debug_export!(PyNcCurrent);

#[pymethods]
impl PyCurrent {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
