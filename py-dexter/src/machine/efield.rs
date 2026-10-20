//! Defines the PyEfield enum that holds one of the Efield objects.

use std::sync::Arc;

use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNoEfield(Arc<NoEfield>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyGaussianEfield(Arc<GaussianEfield>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyEfield", frozen, immutable_type)]
pub enum PyEfield {
    None(PyNoEfield),
    Gaussian(PyGaussianEfield),
}

#[pymethods] // Builders
impl PyEfield {
    #[classmethod]
    pub fn build_no_efield<'py>(_: Bound<'py, PyType>) -> Result<Self> {
        Ok(Self::None(PyNoEfield(Arc::new(NoEfield))))
    }

    #[classmethod]
    pub fn build_gaussian_efield<'py>(
        _: Bound<'py, PyType>,
        strength: f64,
        mean: &PyMagneticFlux,
        width: f64,
    ) -> Result<Self> {
        let inner = GaussianEfield::new(strength, mean.0, width);
        Ok(Self::Gaussian(PyGaussianEfield(Arc::new(inner))))
    }
}

/// References to the trait object and variants
impl PyEfield {
    pub fn inner(&self) -> &dyn Efield {
        match self {
            PyEfield::None(efield) => efield.0.as_ref(),
            PyEfield::Gaussian(efield) => efield.0.as_ref(),
        }
    }

    pub fn no(&self) -> Result<&NoEfield> {
        match self {
            Self::None(efield) => Ok(efield.0.as_ref()),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Efield".into(),
                inner: "NoEfield".into(),
            }),
        }
    }

    pub fn gaussian(&self) -> Result<&GaussianEfield> {
        match self {
            Self::Gaussian(efield) => Ok(efield.0.as_ref()),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Efield".into(),
                inner: "GaussianEfield".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // MachineObject Trait
impl PyEfield {
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

#[pymethods] // Efield Trait
impl PyEfield {
    pub fn eval_e(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_e(flux)?)
    }

    pub fn eval_phi(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_phi(flux)?)
    }

    pub fn eval_phi_deriv(&self, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_phi_deriv(flux)?)
    }
}

// ===============================================================================================

// #[pymethods] // NoEfield
impl PyEfield {}

#[pymethods] // GaussianEfield
impl PyEfield {
    #[getter]
    pub fn strength(&self) -> Result<f64> {
        Ok(self.gaussian()?.strength())
    }

    #[getter]
    pub fn mean(&self) -> Result<PyMagneticFlux> {
        let mean = self.gaussian()?.mean();
        Ok(PyMagneticFlux::from(mean))
    }

    #[getter]
    pub fn width(&self) -> Result<f64> {
        Ok(self.gaussian()?.width())
    }
}

// ===============================================================================================

wrapper_debug_export!(PyNoEfield);
wrapper_debug_export!(PyGaussianEfield);

#[pymethods]
impl PyEfield {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
