//! Defines the PyMode enum that holds one of the Mode objects.

use std::sync::Arc;

use numpy::{IntoPyArray, PyArray1};
use pyo3::{prelude::*, types::PyType};

use crate::*;
use dexter::dexter_machine::*;

// ===============================================================================================

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyFluteMode(Arc<FluteMode>);

#[pyclass(frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub struct PyNcFluteMode(Arc<NcFluteMode>);

// ===============================================================================================

/// Actual export
#[pyclass(name = "_PyMode", frozen, immutable_type, from_py_object)]
#[derive(Clone)]
pub enum PyMode {
    Flute(PyFluteMode),
    Nc(PyNcFluteMode),
}

#[pymethods] // Builders
impl PyMode {
    #[classmethod]
    pub fn build_flute<'py>(
        _: Bound<'py, PyType>,
        epsilon: f64,
        lcfs: &PyMagneticFlux,
        m: i64,
        n: i64,
        phase: f64,
    ) -> Result<Self> {
        let inner = PyFluteMode(Arc::new(FluteMode::new(epsilon, lcfs.0, m, n, phase)));
        Ok(Self::Flute(inner))
    }

    #[classmethod]
    pub fn build_nc<'py>(
        _: Bound<'py, PyType>,
        path: String,
        interp_type: String,
        m: i64,
        n: i64,
        phase_method: Bound<'py, PyAny>,
        analytical_threshold_index: usize,
    ) -> Result<Self> {
        let path = std::path::PathBuf::from(path);
        let typ = resolve_interpolation_1d_type(interp_type)?;
        let phase_method = resolve_phase_method(phase_method)?;
        let builder = NcFluteModeBuilder::new(&path, typ, m, n)
            .with_phase_method(phase_method)
            .with_analytical_threshold_index(analytical_threshold_index);
        let mode = builder.build()?;
        let inner = PyNcFluteMode(Arc::new(mode));

        assert_eq!(Arc::strong_count(&inner.0), 1);
        Ok(Self::Nc(inner))
    }
}

/// References to the trait object and variants
impl PyMode {
    pub fn inner(&self) -> &dyn Mode {
        match self {
            PyMode::Flute(mode) => mode.0.as_ref(),
            PyMode::Nc(mode) => mode.0.as_ref(),
        }
    }

    pub fn boxed_mode(&self) -> Box<dyn Mode> {
        match self {
            PyMode::Flute(mode) => Box::new(Arc::unwrap_or_clone(mode.0.clone())),
            PyMode::Nc(mode) => Box::new(Arc::unwrap_or_clone(mode.0.clone())),
        }
    }

    pub fn flute(&self) -> Result<&FluteMode> {
        match self {
            Self::Flute(mode) => Ok(&mode.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Mode".into(),
                inner: "FluteMode".into(),
            }),
        }
    }

    pub fn nc(&self) -> Result<&NcFluteMode> {
        match self {
            Self::Nc(mode) => Ok(&mode.0),
            _ => Err(DexterError::InvalidVariant {
                wrapper: "Mode".into(),
                inner: "NcFluteMode".into(),
            }),
        }
    }
}

// ===============================================================================================

#[pymethods] // MachineObject Trait
impl PyMode {
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

#[pymethods] // Mode Trait
impl PyMode {
    #[getter]
    pub fn m(&self) -> Result<i64> {
        Ok(self.inner().m())
    }

    #[getter]
    pub fn n(&self) -> Result<i64> {
        Ok(self.inner().n())
    }

    pub fn eval_amplitude(
        &self,
        theta: f64,
        zeta: f64,
        t: f64,
        psi: f64,
        psip: f64,
    ) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self
            .inner()
            .eval_amplitude(flux, theta, zeta, t, &mut self.inner().generate_cache())?)
    }

    pub fn eval_phase(&self, theta: f64, zeta: f64, t: f64, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self
            .inner()
            .eval_phase(flux, theta, zeta, t, &mut self.inner().generate_cache())?)
    }

    pub fn eval_m(&self, theta: f64, zeta: f64, t: f64, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self
            .inner()
            .eval_m(flux, theta, zeta, t, &mut self.inner().generate_cache())?)
    }

    pub fn eval_deriv_flux(
        &self,
        theta: f64,
        zeta: f64,
        t: f64,
        psi: f64,
        psip: f64,
    ) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_deriv_flux(
            flux,
            theta,
            zeta,
            t,
            &mut self.inner().generate_cache(),
        )?)
    }

    pub fn eval_deriv_theta(
        &self,
        theta: f64,
        zeta: f64,
        t: f64,
        psi: f64,
        psip: f64,
    ) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_deriv_theta(
            flux,
            theta,
            zeta,
            t,
            &mut self.inner().generate_cache(),
        )?)
    }

    pub fn eval_deriv_zeta(
        &self,
        theta: f64,
        zeta: f64,
        t: f64,
        psi: f64,
        psip: f64,
    ) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self.inner().eval_deriv_zeta(
            flux,
            theta,
            zeta,
            t,
            &mut self.inner().generate_cache(),
        )?)
    }

    pub fn eval_deriv_t(&self, theta: f64, zeta: f64, t: f64, psi: f64, psip: f64) -> Result<f64> {
        let flux = flux_from_params(psi, psip);
        Ok(self
            .inner()
            .eval_deriv_t(flux, theta, zeta, t, &mut self.inner().generate_cache())?)
    }
}

// ===============================================================================================

#[pymethods] // Flute
impl PyMode {
    #[getter]
    pub fn lcfs(&self) -> Result<PyMagneticFlux> {
        Ok(self.flute()?.lcfs().into())
    }

    #[getter]
    pub fn epsilon(&self) -> Result<f64> {
        Ok(self.flute()?.epsilon())
    }

    #[getter]
    pub fn phase(&self) -> Result<f64> {
        Ok(self.flute()?.phase())
    }
}

#[pymethods] // NcFlute
impl PyMode {
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
    pub fn phase_method(&self) -> Result<String> {
        Ok(format!("{:?}", self.nc()?.phase_method()))
    }

    #[getter]
    pub fn analytical_threshold_index(&self) -> Result<usize> {
        Ok(self.nc()?.analytical_threshold_index())
    }

    #[getter]
    pub fn phase_average(&self) -> Result<Option<f64>> {
        Ok(self.nc()?.phase_average())
    }

    pub fn get_array<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyArray1<f64>>>> {
        let mode = self.nc()?;
        match name {
            "alpha_array" => return Ok(Some(mode.alpha_array().into_pyarray(py))),
            "phase_array" => return Ok(Some(mode.phase_array().into_pyarray(py))),
            "psi_array" => match mode.psi_array() {
                Some(array) => return Ok(Some(array.into_pyarray(py))),
                None => return Ok(None),
            },
            "psip_array" => match mode.psip_array() {
                Some(array) => return Ok(Some(array.into_pyarray(py))),
                None => return Ok(None),
            },
            _ => Err(DexterError::AttributeError {
                obj: "NcBfield".into(),
                attr: name.into(),
            }),
        }
    }
}

// ===============================================================================================

wrapper_debug_export!(PyFluteMode);
wrapper_debug_export!(PyNcFluteMode);

#[pymethods]
impl PyMode {
    pub fn __repr__(&self) -> String {
        format!("{:#?}", self.inner())
    }
}
