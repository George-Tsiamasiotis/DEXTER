//! Defines the PyGeometry enum that holds one of the Geometry objects.

use std::sync::Arc;

use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyLarGeometry(Arc<LarGeometry>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcGeometry(Arc<NcGeometry>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyGeometry", frozen, immutable_type)]
pub enum PyGeometry {
    Lar(PyLarGeometry),
    Nc(PyNcGeometry),
}

#[pymethods] // Builders
impl PyGeometry {
    #[classmethod]
    pub fn build_lar<'py>(
        _: Bound<'py, PyType>,
        baxis: f64,
        raxis: f64,
        rlast: f64,
    ) -> Result<Self> {
        let inner = PyLarGeometry(Arc::new(LarGeometry::new(baxis, raxis, rlast)));
        Ok(Self::Lar(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp1d_type: String,
        interp2d_type: String,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ1 = resolve_interpolation_1d_type(interp1d_type)?;
        let typ2 = resolve_interpolation_2d_type(interp2d_type)?;
        let builder = NcGeometryBuilder::new(&path, typ1, typ2);
        let geometry = builder.build()?;
        let inner = PyNcGeometry(Arc::new(geometry));
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyGeometry {
    pub fn inner(&self) -> &dyn Geometry {
        match self {
            PyGeometry::Lar(geometry) => geometry.0.as_ref(),
            PyGeometry::Nc(geometry) => geometry.0.as_ref(),
        }
    }

    pub fn lar(&self) -> Result<&LarGeometry> {
        match self {
            Self::Lar(geometry) => Ok(&geometry.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Geometry".into(),
                inner: "LarGeometry".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcGeometry> {
        match self {
            Self::Nc(geometry) => Ok(&geometry.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Geometry".into(),
                inner: "NcGeometry".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // MachineObject Trait
impl PyGeometry {
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

#[pymethods] // Geometry Trait
#[rustfmt::skip]
impl PyGeometry {
    #[getter]
    pub fn baxis(&self) -> f64 {
        self.inner().baxis()
    }

    #[getter]
    pub fn raxis(&self) -> f64 {
        self.inner().raxis()
    }

    #[getter]
    pub fn zaxis(&self) -> f64 {
        self.inner().zaxis()
    }

    #[getter]
    pub fn rgeo(&self) -> f64 {
        self.inner().rgeo()
    }

    #[getter]
    pub fn rlast(&self) -> f64 {
        self.inner().rlast()
    }

    pub fn eval_r(&self, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.inner().eval_r(flux, &mut Accelerator::new())?)
    }

    pub fn eval_psit_of_r(&self, r: f64) -> Result<f64> {
        Ok(self.inner().eval_psit_of_r(r, &mut Accelerator::new())?.value())
    }

    pub fn eval_psip_of_r(&self, r: f64) -> Result<f64> {
        Ok(self.inner().eval_psip_of_r(r, &mut Accelerator::new())?.value())
    }

    pub fn eval_rlab(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.inner().eval_rlab(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_zlab(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.inner().eval_zlab(flux, theta, &mut Accelerator2d::new())?)
    }

    pub fn eval_jacobian(&self, theta: f64, psit: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psit, psip);
        Ok(self.inner().eval_jacobian(flux, theta, &mut Accelerator2d::new())?)
    }

    #[getter]
    pub fn rlab_last<'py>(&self, py: Python<'py>)  -> Result<Bound<'py, PyArray1<f64>>> {
        Ok(self.inner().rlab_last().into_pyarray(py))
    }

    #[getter]
    pub fn zlab_last<'py>(&self, py: Python<'py>)  -> Result<Bound<'py, PyArray1<f64>>> {
        Ok(self.inner().zlab_last().into_pyarray(py))
    }
}

// ===============================================================================================

#[pymethods] // "Common" but not provided by the `Geometry` trait
impl PyGeometry {
    #[getter]
    pub fn psit_last(&self) -> Option<PyMagneticFlux> {
        match self {
            Self::Lar(inner) => Some(inner.0.psit_last().into()), // guaranteed to exist
            Self::Nc(inner) => inner.0.psit_last().map(Into::into),
        }
    }

    #[getter]
    pub fn psip_last(&self) -> Option<PyMagneticFlux> {
        match self {
            Self::Lar(_) => None, // guaranteed to **not** exist
            Self::Nc(inner) => inner.0.psip_last().map(Into::into),
        }
    }
}

#[pymethods] // Lar
impl PyGeometry {}

#[pymethods] // Nc
impl PyGeometry {
    #[getter]
    pub fn path(&self) -> Result<String> {
        Ok(self.nc()?.path().to_str().unwrap_or_default().to_string())
    }

    #[getter]
    pub fn netcdf_version(&self) -> Result<String> {
        Ok(self.nc()?.netcdf_version().to_string())
    }

    #[getter]
    pub fn interp1d_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp1d_type()))
    }

    #[getter]
    pub fn interp2d_type(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.interp2d_type()))
    }

    #[getter]
    pub fn shape(&self) -> Result<(usize, usize)> {
        Ok(self.nc()?.shape())
    }

    pub fn get_array<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        let geometry = self.nc()?;
        let values: &[f64] = match name {
            "theta_array" => geometry.theta_values(),
            "r_array" => geometry.r_values(),
            "psit_array" => match geometry.psit_values() {
                Some(values) => values,
                None => return Ok(None),
            },
            "psip_array" => match geometry.psip_values() {
                Some(values) => values,
                None => return Ok(None),
            },
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "NcGeometry".into(),
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
    ) -> Result<Option<Bound<'py, PyArray2<f64>>>> {
        let geometry = self.nc()?;
        let view = match name {
            "rlab_array" => geometry.rlab_array(),
            "zlab_array" => geometry.zlab_array(),
            "jacobian_array" => match geometry.jacobian_array() {
                Some(array) => array,
                None => return Ok(None),
            },
            _ => {
                return Err(DexterError::AttributeError {
                    obj: "NcGeometry".into(),
                    attr: name.into(),
                });
            }
        };
        Ok(Some(view.to_owned().into_pyarray(py)))
    }
}

// ===============================================================================================

wrapper_debug_export!(PyLarGeometry);
wrapper_debug_export!(PyNcGeometry);

#[pymethods]
impl PyGeometry {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
