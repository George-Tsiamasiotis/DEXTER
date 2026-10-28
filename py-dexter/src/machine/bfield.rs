//! Defines the PyBfield enum that holds one of the Bfield objects.

use std::sync::Arc;

use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyLarBfield(Arc<LarBfield>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcBfield(Arc<NcBfield>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyBfield", frozen, immutable_type)]
pub enum PyBfield {
    Lar(PyLarBfield),
    Nc(PyNcBfield),
}

#[pymethods] // Builders
impl PyBfield {
    #[classmethod]
    pub fn build_lar<'py>(_: Bound<'py, PyType>) -> Result<Self> {
        let inner = PyLarBfield(Arc::new(LarBfield));
        Ok(Self::Lar(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp_type: String,
        padding: usize,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_2d_type(interp_type)?;
        let builder = NcBfieldBuilder::new(&path, typ).with_padding(padding);
        let bfield = builder.build()?;
        let inner = PyNcBfield(Arc::new(bfield));
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyBfield {
    pub fn inner(&self) -> &dyn Bfield {
        match self {
            PyBfield::Lar(bfield) => bfield.0.as_ref(),
            PyBfield::Nc(bfield) => bfield.0.as_ref(),
        }
    }

    pub fn lar(&self) -> Result<&LarBfield> {
        match self {
            Self::Lar(bfield) => Ok(&bfield.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Bfield".into(),
                inner: "LarBfield".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcBfield> {
        match self {
            Self::Nc(bfield) => Ok(&bfield.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Bfield".into(),
                inner: "NcBfield".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // MachineObject Trait
impl PyBfield {
    #[getter]
    pub fn machine_type(&self) -> String {
        format!("{:?}", self.inner().machine_type())
    }

    #[getter]
    pub fn psit_state(&self) -> String {
        format!("{:?}", self.inner().psit_state())
    }

    #[getter]
    pub fn psip_state(&self) -> String {
        format!("{:?}", self.inner().psip_state())
    }
}

#[pymethods] // Bfield Trait
impl PyBfield {
    pub fn eval_b(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_b(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_deriv_flux(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_deriv_flux(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_deriv_theta(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_deriv_theta(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_deriv_flux2(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_deriv_flux2(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_deriv_theta2(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_deriv_theta2(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_deriv_mixed(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self
            .inner()
            .eval_deriv_mixed(flux, theta, &mut Accelerator2d::new())?)
    }
}

// ===============================================================================================

#[pymethods] // Lar
impl PyBfield {}

#[pymethods] // Nc
impl PyBfield {
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

    #[getter]
    pub fn baxis(&self) -> Result<f64> {
        Ok(self.nc()?.baxis())
    }

    #[getter]
    pub fn padding(&self) -> Result<usize> {
        Ok(self.nc()?.padding())
    }

    #[getter]
    pub fn padding_theta(&self) -> Result<f64> {
        Ok(self.nc()?.padding_theta())
    }

    #[getter]
    pub fn shape(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape())
    }

    #[getter]
    pub fn shape_padded(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape_padded())
    }

    pub fn get_array<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        let bfield = self.nc()?;
        let values: &[f64] = match name {
            "theta_array" => bfield.theta_values(),
            "theta_array_padded" => bfield.theta_array_padded(),
            "psit_array" => match bfield.psit_values() {
                Some(values) => values,
                None => return Ok(None),
            },
            "psip_array" => match bfield.psip_values() {
                Some(values) => values,
                None => return Ok(None),
            },
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "NcBfield".into(),
                    attr: name.into(),
                });
            }
        };
        Ok(Some(Array1::from(values.to_vec()).into_pyarray(py)))
    }

    pub fn get_array2d<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Bound<'py, PyArray2<f64>>> {
        let bfield = self.nc()?;
        let view = match name {
            "b_array" => bfield.b_array(),
            "b_array_padded" => bfield.b_array_padded(),
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "NcBfield".into(),
                    attr: name.into(),
                });
            }
        };
        Ok(view.to_owned().into_pyarray(py))
    }
}

// ===============================================================================================

wrapper_debug_export!(PyLarBfield);
wrapper_debug_export!(PyNcBfield);

#[pymethods]
impl PyBfield {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
