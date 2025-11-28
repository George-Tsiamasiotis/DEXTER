//! Equilibrium objects' Python wrappers.

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyList;
use rsl_interpolation::{Accelerator, Cache};
use safe_unwrap::safe_unwrap;

use equilibrium::{Bfield, Currents, Geometry, Harmonic, HarmonicCache, Perturbation, Qfactor};
use equilibrium::{Flux, Length};
use utils::*;

use super::pyerrors::PyEqError;

// ===============================================================================================

#[pyclass(frozen, name = "Geometry")]
pub struct PyGeometry(pub Geometry);

#[pymethods]
impl PyGeometry {
    #[new]
    pub fn new(path: &str, typ1d: &str, typ2d: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Geometry::from_dataset(&path, typ1d, typ2d)?))
    }

    #[getter]
    pub fn get_shape(&self) -> (usize, usize) {
        self.0.shape()
    }
}

py_debug_impl!(PyGeometry);
py_repr_impl!(PyGeometry);
py_get_path!(PyGeometry);
py_export_getter!(PyGeometry, typ1d, String);
py_export_getter!(PyGeometry, typ2d, String);
py_export_getter!(PyGeometry, baxis, f64);
py_export_getter!(PyGeometry, raxis, f64);
py_export_getter!(PyGeometry, psip_wall, Flux);
py_export_getter!(PyGeometry, psi_wall, Flux);
py_export_getter!(PyGeometry, r_wall, Length);
py_get_numpy1D!(PyGeometry, theta_data);
py_get_numpy1D!(PyGeometry, psip_data);
py_get_numpy1D!(PyGeometry, psi_data);
py_get_numpy1D!(PyGeometry, r_data);
py_get_numpy2D!(PyGeometry, rlab_data);
py_get_numpy2D!(PyGeometry, zlab_data);
py_eval1D!(PyGeometry, r, "no-acc");
py_eval1D!(PyGeometry, psip, "no-acc");
py_eval2D!(PyGeometry, rlab, "no-acc");
py_eval2D!(PyGeometry, zlab, "no-acc");

// ===============================================================================================

#[pyclass(frozen, name = "Qfactor")]
pub struct PyQfactor(pub Qfactor);

#[pymethods]
impl PyQfactor {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Qfactor::from_dataset(&path, typ)?))
    }
}

py_debug_impl!(PyQfactor);
py_repr_impl!(PyQfactor);
py_get_path!(PyQfactor);
py_len_impl!(PyQfactor);
py_export_getter!(PyQfactor, typ, String);
py_export_getter!(PyQfactor, psip_wall, f64);
py_export_getter!(PyQfactor, psi_wall, f64);
py_get_numpy1D!(PyQfactor, psip_data);
py_get_numpy1D!(PyQfactor, q_data);
py_get_numpy1D!(PyQfactor, r_data);
py_get_numpy1D!(PyQfactor, psi_data);
py_get_numpy1D_fallible!(PyQfactor, q_data_derived);
py_eval1D!(PyQfactor, q);
py_eval1D!(PyQfactor, psi);

// ===============================================================================================

#[pyclass(frozen, name = "Currents")]
pub struct PyCurrents(pub Currents);

#[pymethods]
impl PyCurrents {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Currents::from_dataset(&path, typ)?))
    }
}

py_debug_impl!(PyCurrents);
py_repr_impl!(PyCurrents);
py_get_path!(PyCurrents);
py_len_impl!(PyCurrents);
py_export_getter!(PyCurrents, typ, String);
py_export_getter!(PyCurrents, psip_wall, f64);
py_eval1D!(PyCurrents, g);
py_eval1D!(PyCurrents, i);
py_eval1D!(PyCurrents, dg_dpsip);
py_eval1D!(PyCurrents, di_dpsip);
py_get_numpy1D!(PyCurrents, psip_data);
py_get_numpy1D!(PyCurrents, g_data);
py_get_numpy1D!(PyCurrents, i_data);

// ===============================================================================================

#[pyclass(frozen, name = "Bfield")]
pub struct PyBfield(pub Bfield);

#[pymethods]
impl PyBfield {
    #[new]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Bfield::from_dataset(&path, typ)?))
    }

    #[getter]
    pub fn get_shape(&self) -> (usize, usize) {
        self.0.shape()
    }
}

py_debug_impl!(PyBfield);
py_repr_impl!(PyBfield);
py_get_path!(PyBfield);
py_export_getter!(PyBfield, typ, String);
py_export_getter!(PyBfield, psip_wall, f64);
py_eval2D!(PyBfield, b);
py_eval2D!(PyBfield, db_dtheta);
py_eval2D!(PyBfield, db_dpsip);
py_eval2D!(PyBfield, d2b_dtheta2);
py_eval2D!(PyBfield, d2b_dpsip2);
py_eval2D!(PyBfield, d2b_dpsip_dtheta);
py_get_numpy1D!(PyBfield, psip_data);
py_get_numpy1D!(PyBfield, theta_data);
py_get_numpy2D!(PyBfield, b_data);
py_get_numpy2D!(PyBfield, db_dpsip_data);
py_get_numpy2D!(PyBfield, db_dtheta_data);

// ===============================================================================================

#[derive(Clone)]
#[pyclass(frozen, name = "Harmonic")]
pub struct PyHarmonic(pub Harmonic);

#[pymethods]
impl PyHarmonic {
    #[new]
    pub fn new(path: &str, typ: &str, m: i64, n: i64) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        Ok(Self(Harmonic::from_dataset(&path, typ, m, n)?))
    }
}

impl From<&Harmonic> for PyHarmonic {
    fn from(harmonic: &Harmonic) -> Self {
        safe_unwrap!(
            "If harmonic exists, Pyharmonic::new() cannot fail",
            PyHarmonic::new(
                safe_unwrap!("file already opened", harmonic.path().to_str()),
                harmonic.typ().as_str(),
                harmonic.m(),
                harmonic.n(),
            )
        )
    }
}

py_debug_impl!(PyHarmonic);
py_repr_impl!(PyHarmonic);
py_get_path!(PyHarmonic);
py_len_impl!(PyHarmonic);
py_export_getter!(PyHarmonic, typ, String);
py_export_getter!(PyHarmonic, psip_wall, f64);
py_export_getter!(PyHarmonic, phase_average, f64);
py_export_getter!(PyHarmonic, m, i64);
py_export_getter!(PyHarmonic, n, i64);
py_get_numpy1D!(PyHarmonic, psip_data);
py_get_numpy1D!(PyHarmonic, a_data);
py_get_numpy1D!(PyHarmonic, phase_data);
py_eval_harmonic!(PyHarmonic, h);
py_eval_harmonic!(PyHarmonic, dh_dpsip);
py_eval_harmonic!(PyHarmonic, dh_dtheta);
py_eval_harmonic!(PyHarmonic, dh_dzeta);
py_eval_harmonic!(PyHarmonic, dh_dt);
py_eval1D!(PyHarmonic, a);
py_eval1D!(PyHarmonic, da_dpsip);
py_eval1D!(PyHarmonic, phase);

// ===============================================================================================

#[pyclass(frozen, name = "Perturbation")]
pub struct PyPerturbation(pub Perturbation);

#[pymethods]
impl PyPerturbation {
    #[new]
    pub fn new_py<'py>(harmonics: Bound<'py, PyList>) -> Result<Self, PyEqError> {
        let pyharmonics_vec: Vec<PyHarmonic> = harmonics
            .iter()
            .map(|ph| {
                ph.extract()
                    .expect("Could not extract 'PyHarmonic' from python list")
            })
            .collect();
        let harmonics_vec: Vec<Harmonic> = pyharmonics_vec
            .clone()
            .iter()
            .map(|ph| ph.0.clone())
            .collect();

        Ok(Self(Perturbation::from_harmonics(&harmonics_vec)))
    }

    /// Makes PyPerturbation indexable
    pub fn __getitem__(&self, index: usize) -> PyHarmonic {
        PyHarmonic::from(
            self.0
                .get_harmonics()
                .get(index)
                .expect("Harmonic index out of bounds"),
        )
    }
}

py_debug_impl!(PyPerturbation);
py_repr_impl!(PyPerturbation);
py_len_impl!(PyPerturbation);
py_eval_perturbation!(PyPerturbation, p);
py_eval_perturbation!(PyPerturbation, dp_dpsip);
py_eval_perturbation!(PyPerturbation, dp_dtheta);
py_eval_perturbation!(PyPerturbation, dp_dzeta);
py_eval_perturbation!(PyPerturbation, dp_dt);
