//! Equilibrium objects' Python wrappers.

use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyList;
use rsl_interpolation::{Accelerator, Cache};

use dexter::equilibrium::*;
use std::result::Result; // awful; replace equilibrium's Result

use super::pyerrors::PyEqError;
use crate::{
    py_debug_impl, py_eval_harmonic, py_eval_perturbation, py_eval1D, py_eval2D, py_export_getter,
    py_get_numpy1D, py_get_numpy2D, py_get_path, py_len_impl, py_repr_impl,
};

// ===============================================================================================
// ===============================================================================================

#[pyclass(frozen, name = "NcGeometry")]
pub struct PyNcGeometry(pub NcGeometry);

#[pymethods]
impl PyNcGeometry {
    #[new]
    #[pyo3(signature = (path, typ1d, typ2d))]
    pub fn new(path: &str, typ1d: &str, typ2d: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
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
py_export_getter!(PyNcGeometry, typ1d, String);
py_export_getter!(PyNcGeometry, typ2d, String);
py_export_getter!(PyNcGeometry, baxis, f64);
py_export_getter!(PyNcGeometry, raxis, f64);
py_export_getter!(PyNcGeometry, zaxis, f64);
py_export_getter!(PyNcGeometry, rgeo, f64);
py_export_getter!(PyNcGeometry, psip_wall, f64);
py_export_getter!(PyNcGeometry, psi_wall, f64);
py_export_getter!(PyNcGeometry, r_wall, f64);
py_get_numpy1D!(PyNcGeometry, theta_data);
py_get_numpy1D!(PyNcGeometry, psip_data);
py_get_numpy1D!(PyNcGeometry, psi_data);
py_get_numpy1D!(PyNcGeometry, r_data);
py_get_numpy2D!(PyNcGeometry, rlab_data);
py_get_numpy2D!(PyNcGeometry, zlab_data);
py_get_numpy2D!(PyNcGeometry, jacobian_data);
py_eval1D!(PyNcGeometry, r, "no-acc");
py_eval1D!(PyNcGeometry, psi, "no-acc");
py_eval1D!(PyNcGeometry, psip, "no-acc");
py_eval2D!(PyNcGeometry, rlab, "no-acc");
py_eval2D!(PyNcGeometry, zlab, "no-acc");
py_eval2D!(PyNcGeometry, jacobian, "no-acc");

// ===============================================================================================
// ===============================================================================================

#[pyclass(frozen, name = "NcQfactor")]
pub struct PyNcQfactor(pub NcQfactor);

#[pymethods]
impl PyNcQfactor {
    #[new]
    #[pyo3(signature = (path, typ))]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcQfactorBuilder::new(&path, typ);
        Ok(Self(builder.build()?))
    }
}

py_debug_impl!(PyNcQfactor);
py_repr_impl!(PyNcQfactor);
py_get_path!(PyNcQfactor);
py_len_impl!(PyNcQfactor);
py_export_getter!(PyNcQfactor, typ, String);
py_get_numpy1D!(PyNcQfactor, psip_data);
py_get_numpy1D!(PyNcQfactor, q_data);
py_get_numpy1D!(PyNcQfactor, psi_data);
py_eval1D!(PyNcQfactor, q);
py_eval1D!(PyNcQfactor, psi);
py_eval1D!(PyNcQfactor, dpsi_dpsip);

// ===============================================================================================

#[pyclass(frozen, name = "UnityQfactor")]
pub struct PyUnityQfactor(pub UnityQfactor);

#[pymethods]
impl PyUnityQfactor {
    #[new]
    #[pyo3(signature = ())]
    pub fn new() -> Result<Self, PyEqError> {
        Ok(Self(UnityQfactor))
    }
}

py_debug_impl!(PyUnityQfactor);
py_repr_impl!(PyUnityQfactor);
py_eval1D!(PyUnityQfactor, q);
py_eval1D!(PyUnityQfactor, psi);
py_eval1D!(PyUnityQfactor, dpsi_dpsip);

// ===============================================================================================
// ===============================================================================================

#[pyclass(frozen, name = "NcCurrent")]
pub struct PyNcCurrent(pub NcCurrent);

#[pymethods]
impl PyNcCurrent {
    #[new]
    #[pyo3(signature = (path, typ))]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcCurrentBuilder::new(&path, typ);
        Ok(Self(builder.build()?))
    }
}

py_debug_impl!(PyNcCurrent);
py_repr_impl!(PyNcCurrent);
py_get_path!(PyNcCurrent);
py_len_impl!(PyNcCurrent);
py_export_getter!(PyNcCurrent, typ, String);
py_eval1D!(PyNcCurrent, g);
py_eval1D!(PyNcCurrent, i);
py_eval1D!(PyNcCurrent, dg_dpsip);
py_eval1D!(PyNcCurrent, di_dpsip);
py_get_numpy1D!(PyNcCurrent, psip_data);
py_get_numpy1D!(PyNcCurrent, g_data);
py_get_numpy1D!(PyNcCurrent, i_data);

// ===============================================================================================

#[pyclass(frozen, name = "LarCurrent")]
pub struct PyLarCurrent(pub LarCurrent);

#[pymethods]
impl PyLarCurrent {
    #[new]
    #[pyo3(signature = ())]
    pub fn new() -> Result<Self, PyEqError> {
        Ok(Self(LarCurrent))
    }
}

py_debug_impl!(PyLarCurrent);
py_repr_impl!(PyLarCurrent);
py_eval1D!(PyLarCurrent, g);
py_eval1D!(PyLarCurrent, i);
py_eval1D!(PyLarCurrent, dg_dpsip);
py_eval1D!(PyLarCurrent, di_dpsip);

// ===============================================================================================

#[pyclass(frozen, name = "NcBfield")]
pub struct PyNcBfield(pub NcBfield);

#[pymethods]
impl PyNcBfield {
    #[new]
    #[pyo3(signature = (path, typ))]
    pub fn new(path: &str, typ: &str) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder = NcBfieldBuilder::new(&path, typ);
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
py_export_getter!(PyNcBfield, typ, String);
py_eval2D!(PyNcBfield, b);
py_eval2D!(PyNcBfield, db_dtheta);
py_eval2D!(PyNcBfield, db_dpsip);
py_get_numpy1D!(PyNcBfield, psip_data);
py_get_numpy1D!(PyNcBfield, theta_data);
py_get_numpy2D!(PyNcBfield, b_data);
py_get_numpy2D!(PyNcBfield, db_dpsip_data);
py_get_numpy2D!(PyNcBfield, db_dtheta_data);

// ===============================================================================================

#[derive(Clone)]
#[pyclass(frozen, name = "NcHarmonic")]
pub struct PyNcHarmonic(pub NcHarmonic);

#[pymethods]
impl PyNcHarmonic {
    #[new]
    #[pyo3(signature = (path, typ, m, n, phase_method = "Resonance"))]
    pub fn new(
        path: &str,
        typ: &str,
        m: i64,
        n: i64,
        phase_method: Option<&str>,
    ) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let builder =
            NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(match phase_method {
                Some(method) => match method.to_lowercase().as_str() {
                    "zero" => PhaseMethod::Zero,
                    "average" => PhaseMethod::Average,
                    "resonance" => PhaseMethod::Resonance,
                    "interpolation" => PhaseMethod::Interpolation,
                    "custom" => todo!("How to pass this?"), // TODO:
                    _ => panic!("Invalid phase method"),
                },
                None => PhaseMethod::default(),
            });

        Ok(Self(builder.build()?))
    }

    #[getter]
    fn get_phase_method(&self) -> String {
        format!("{:?}", self.0.phase_method())
    }
}

impl From<&NcHarmonic> for PyNcHarmonic {
    fn from(harmonic: &NcHarmonic) -> Self {
        PyNcHarmonic::new(
            harmonic.path().to_str().unwrap(), // Safe: already exists
            harmonic.typ().as_str(),
            harmonic.m(),
            harmonic.n(),
            Some(match harmonic.phase_method() {
                PhaseMethod::Zero => "zero",
                PhaseMethod::Average => "average",
                PhaseMethod::Resonance => "resonance",
                PhaseMethod::Interpolation => "interpolation",
                PhaseMethod::Custom(_) => "custom",
            }),
        )
        .unwrap()
    }
}

py_debug_impl!(PyNcHarmonic);
py_repr_impl!(PyNcHarmonic);
py_get_path!(PyNcHarmonic);
py_len_impl!(PyNcHarmonic);
py_export_getter!(PyNcHarmonic, typ, String);
py_export_getter!(PyNcHarmonic, phase_average, Option<f64>);
py_export_getter!(PyNcHarmonic, phase_resonance, Option<f64>);
py_export_getter!(PyNcHarmonic, m, i64);
py_export_getter!(PyNcHarmonic, n, i64);
py_get_numpy1D!(PyNcHarmonic, psip_data);
py_get_numpy1D!(PyNcHarmonic, a_data);
py_get_numpy1D!(PyNcHarmonic, phase_data);
py_eval_harmonic!(PyNcHarmonic, h);
py_eval_harmonic!(PyNcHarmonic, dh_dpsip);
py_eval_harmonic!(PyNcHarmonic, dh_dtheta);
py_eval_harmonic!(PyNcHarmonic, dh_dzeta);
py_eval_harmonic!(PyNcHarmonic, dh_dt);
py_eval1D!(PyNcHarmonic, a);
py_eval1D!(PyNcHarmonic, da_dpsip);
py_eval1D!(PyNcHarmonic, phase);

// ===============================================================================================

#[pyclass(frozen, name = "NcPerturbation")]
pub struct PyNcPerturbation(pub NcPerturbation);

#[pymethods]
impl PyNcPerturbation {
    #[new]
    #[pyo3(signature = (harmonics))]
    pub fn new_py<'py>(harmonics: Bound<'py, PyList>) -> Result<Self, PyEqError> {
        let pyharmonics_vec: Vec<PyNcHarmonic> = harmonics
            .iter()
            .map(|ph| {
                ph.extract()
                    .expect("Could not extract 'PyNcHarmonic' from python list")
            })
            .collect();
        let harmonics_vec: Vec<NcHarmonic> = pyharmonics_vec
            .clone()
            .iter()
            .map(|ph| ph.0.clone())
            .collect();

        Ok(Self(NcPerturbation::from_harmonics(&harmonics_vec)))
    }

    /// Makes PyNcPerturbation indexable
    pub fn __getitem__(&self, index: usize) -> PyNcHarmonic {
        PyNcHarmonic::from(
            self.0
                .get_harmonics()
                .get(index)
                .expect("NcHarmonic index out of bounds"),
        )
    }
}

py_debug_impl!(PyNcPerturbation);
py_repr_impl!(PyNcPerturbation);
py_len_impl!(PyNcPerturbation);
py_eval_perturbation!(PyNcPerturbation, p);
py_eval_perturbation!(PyNcPerturbation, dp_dpsip);
py_eval_perturbation!(PyNcPerturbation, dp_dtheta);
py_eval_perturbation!(PyNcPerturbation, dp_dzeta);
py_eval_perturbation!(PyNcPerturbation, dp_dt);
