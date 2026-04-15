//! `dexter-equilibrium` newtypes, constructors and method exports.

use dexter::dexter_equilibrium::Perturbation;
use dexter::dexter_equilibrium::{Bfield, LarBfield, NcBfield, NcBfieldBuilder};
use dexter::dexter_equilibrium::{
    CosHarmonic, Harmonic, NcHarmonic, NcHarmonicBuilder, PhaseMethod,
};
use dexter::dexter_equilibrium::{Current, LarCurrent, NcCurrent, NcCurrentBuilder};
use dexter::dexter_equilibrium::{FluxCommute, LastClosedFluxSurface};
use dexter::dexter_equilibrium::{Geometry, LarGeometry, NcGeometry, NcGeometryBuilder};
use dexter::dexter_equilibrium::{
    NcQfactor, NcQfactorBuilder, ParabolicQfactor, Qfactor, UnityQfactor,
};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyInt, PyList, PyTuple};
use rsl_interpolation::{Accelerator, Cache};

use crate::pyerror::{PyEqError, PyEvalError};
use crate::{
    py_debug_impl, py_eval_harmonic, py_eval_perturbation, py_eval1D, py_eval2D, py_export_getter,
    py_get_enum_string, py_get_netcdf_version, py_get_numpy1D, py_get_numpy1D_fallible,
    py_get_numpy2D, py_get_path, py_repr_impl,
};

#[pyclass(name = "_PyLastClosedFluxSurface", frozen)]
pub struct PyLastClosedFluxSurface(pub(crate) LastClosedFluxSurface);

#[pymethods]
impl PyLastClosedFluxSurface {
    #[new]
    pub fn new(kind: &str, value: f64) -> PyResult<Self> {
        Ok(Self(match kind.to_lowercase().as_str() {
            "toroidal" => LastClosedFluxSurface::Toroidal(value),
            "poloidal" => LastClosedFluxSurface::Poloidal(value),
            _ => {
                return Err(PyErr::new::<PyTypeError, _>(
                    "Invalid 'LastClosedFluxSurface'",
                ));
            }
        }))
    }

    #[getter]
    pub fn get_kind(&self) -> &str {
        match self.0 {
            LastClosedFluxSurface::Toroidal(_) => "Toroidal",
            LastClosedFluxSurface::Poloidal(_) => "Poloidal",
        }
    }

    #[getter]
    pub fn get_value(&self) -> f64 {
        self.0.value()
    }
}

py_debug_impl!(PyLastClosedFluxSurface);
py_repr_impl!(PyLastClosedFluxSurface);

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "_PyLarGeometry", frozen)]
pub struct PyLarGeometry(pub(crate) LarGeometry);

#[pymethods]
impl PyLarGeometry {
    #[new]
    #[pyo3(signature = (baxis, raxis, rlast))]
    pub fn new(baxis: f64, raxis: f64, rlast: f64) -> Self {
        Self(LarGeometry::new(baxis, raxis, rlast))
    }
}

py_debug_impl!(PyLarGeometry);
py_repr_impl!(PyLarGeometry);
py_get_enum_string!(PyLarGeometry, equilibrium_type);
py_export_getter!(PyLarGeometry, baxis, f64);
py_export_getter!(PyLarGeometry, raxis, f64);
py_export_getter!(PyLarGeometry, zaxis, f64);
py_export_getter!(PyLarGeometry, rgeo, f64);
py_export_getter!(PyLarGeometry, rlast, f64);
py_export_getter!(PyLarGeometry, psi_last, f64);
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
py_get_numpy1D!(PyLarGeometry, rlab_last);
py_get_numpy1D!(PyLarGeometry, zlab_last);

// ===============================================================================================

#[pyclass(name = "_PyNcGeometry", frozen)]
pub struct PyNcGeometry(pub(crate) NcGeometry);

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
py_export_getter!(PyNcGeometry, rlast, f64);
py_export_getter!(PyNcGeometry, psi_last, Option<f64>);
py_export_getter!(PyNcGeometry, psip_last, Option<f64>);
py_get_enum_string!(PyNcGeometry, psi_state);
py_get_enum_string!(PyNcGeometry, psip_state);
py_get_numpy1D_fallible!(PyNcGeometry, psi_array);
py_get_numpy1D_fallible!(PyNcGeometry, psip_array);
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
py_get_numpy1D!(PyNcGeometry, rlab_last);
py_get_numpy1D!(PyNcGeometry, zlab_last);

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "_PyUnityQfactor", frozen)]
pub struct PyUnityQfactor(pub(crate) UnityQfactor);

#[pymethods]
impl PyUnityQfactor {
    #[new]
    #[pyo3(signature = (lcfs))]
    pub fn new(lcfs: &PyLastClosedFluxSurface) -> Self {
        Self(UnityQfactor::new(lcfs.0))
    }
}

py_debug_impl!(PyUnityQfactor);
py_repr_impl!(PyUnityQfactor);
py_get_enum_string!(PyUnityQfactor, equilibrium_type);
py_export_getter!(PyUnityQfactor, psi_last, f64);
py_export_getter!(PyUnityQfactor, psip_last, f64);
py_export_getter!(PyUnityQfactor, qlast, f64);
py_export_getter!(PyUnityQfactor, qaxis, f64);
py_eval1D!(PyUnityQfactor, psip_of_psi);
py_eval1D!(PyUnityQfactor, psi_of_psip);
py_eval1D!(PyUnityQfactor, q_of_psi);
py_eval1D!(PyUnityQfactor, q_of_psip);
py_eval1D!(PyUnityQfactor, dpsi_dpsip);
py_eval1D!(PyUnityQfactor, dpsip_dpsi);
py_eval1D!(PyUnityQfactor, iota_of_psi);
py_eval1D!(PyUnityQfactor, iota_of_psip);
py_eval1D!(PyUnityQfactor, psi_of_q);
py_eval1D!(PyUnityQfactor, psip_of_q);

// ===============================================================================================

#[pyclass(name = "_PyParabolicQfactor", frozen)]
pub struct PyParabolicQfactor(pub(crate) ParabolicQfactor);

#[pymethods]
impl PyParabolicQfactor {
    #[new]
    #[pyo3(signature = (qaxis, qlast, lcfs))]
    pub fn new<'py>(qaxis: f64, qlast: f64, lcfs: &PyLastClosedFluxSurface) -> PyResult<Self> {
        Ok(Self(ParabolicQfactor::new(qaxis, qlast, lcfs.0)))
    }
}

py_debug_impl!(PyParabolicQfactor);
py_repr_impl!(PyParabolicQfactor);
py_get_enum_string!(PyParabolicQfactor, equilibrium_type);
py_export_getter!(PyParabolicQfactor, psi_last, f64);
py_export_getter!(PyParabolicQfactor, psip_last, f64);
py_export_getter!(PyParabolicQfactor, qlast, f64);
py_export_getter!(PyParabolicQfactor, qaxis, f64);
py_eval1D!(PyParabolicQfactor, psip_of_psi);
py_eval1D!(PyParabolicQfactor, psi_of_psip);
py_eval1D!(PyParabolicQfactor, q_of_psi);
py_eval1D!(PyParabolicQfactor, q_of_psip);
py_eval1D!(PyParabolicQfactor, dpsi_dpsip);
py_eval1D!(PyParabolicQfactor, dpsip_dpsi);
py_eval1D!(PyParabolicQfactor, iota_of_psi);
py_eval1D!(PyParabolicQfactor, iota_of_psip);
py_eval1D!(PyParabolicQfactor, psi_of_q);
py_eval1D!(PyParabolicQfactor, psip_of_q);

// ===============================================================================================

#[pyclass(name = "_PyNcQfactor", frozen)]
pub struct PyNcQfactor(pub(crate) NcQfactor);

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
py_export_getter!(PyNcQfactor, psi_last, f64);
py_export_getter!(PyNcQfactor, psip_last, f64);
py_export_getter!(PyNcQfactor, qlast, f64);
py_export_getter!(PyNcQfactor, qaxis, f64);
py_export_getter!(PyNcQfactor, interp_type, String);
py_get_enum_string!(PyNcQfactor, psi_state);
py_get_enum_string!(PyNcQfactor, psip_state);
py_get_numpy1D_fallible!(PyNcQfactor, psi_array);
py_get_numpy1D_fallible!(PyNcQfactor, psip_array);
py_get_numpy1D!(PyNcQfactor, q_array);
py_eval1D!(PyNcQfactor, psip_of_psi);
py_eval1D!(PyNcQfactor, psi_of_psip);
py_eval1D!(PyNcQfactor, q_of_psi);
py_eval1D!(PyNcQfactor, q_of_psip);
py_eval1D!(PyNcQfactor, dpsi_dpsip);
py_eval1D!(PyNcQfactor, dpsip_dpsi);
py_eval1D!(PyNcQfactor, iota_of_psi);
py_eval1D!(PyNcQfactor, iota_of_psip);
py_eval1D!(PyNcQfactor, psi_of_q);
py_eval1D!(PyNcQfactor, psip_of_q);

// ===============================================================================================
// ===============================================================================================

#[pyclass(name = "_PyLarCurrent", frozen)]
pub struct PyLarCurrent(pub(crate) LarCurrent);

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

#[pyclass(name = "_PyNcCurrent", frozen)]
pub struct PyNcCurrent(pub(crate) NcCurrent);

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
py_export_getter!(PyNcCurrent, psi_last, Option<f64>);
py_export_getter!(PyNcCurrent, psip_last, Option<f64>);
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

// ===============================================================================================
// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyCosHarmonic", frozen)]
pub struct PyCosHarmonic(pub(crate) CosHarmonic);

#[pymethods]
impl PyCosHarmonic {
    #[new]
    #[pyo3(signature = (epsilon, lcfs, m, n, phase))]
    pub fn new(epsilon: f64, lcfs: &PyLastClosedFluxSurface, m: i64, n: i64, phase: f64) -> Self {
        Self(CosHarmonic::new(epsilon, lcfs.0, m, n, phase))
    }

    #[getter]
    pub fn get_lcfs(&self) -> PyLastClosedFluxSurface {
        PyLastClosedFluxSurface(self.0.lcfs())
    }
}

py_debug_impl!(PyCosHarmonic);
py_repr_impl!(PyCosHarmonic);
py_get_enum_string!(PyCosHarmonic, equilibrium_type);
py_export_getter!(PyCosHarmonic, epsilon, f64);
py_export_getter!(PyCosHarmonic, psi_last, Option<f64>);
py_export_getter!(PyCosHarmonic, psip_last, Option<f64>);
py_export_getter!(PyCosHarmonic, phase, f64);
py_export_getter!(PyCosHarmonic, m, i64);
py_export_getter!(PyCosHarmonic, n, i64);
py_eval_harmonic!(PyCosHarmonic, alpha_of_psi);
py_eval_harmonic!(PyCosHarmonic, alpha_of_psip);
py_eval_harmonic!(PyCosHarmonic, phase_of_psi);
py_eval_harmonic!(PyCosHarmonic, phase_of_psip);
py_eval_harmonic!(PyCosHarmonic, h_of_psi);
py_eval_harmonic!(PyCosHarmonic, h_of_psip);
py_eval_harmonic!(PyCosHarmonic, dh_dpsi);
py_eval_harmonic!(PyCosHarmonic, dh_dpsip);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dtheta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dtheta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dzeta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dzeta);
py_eval_harmonic!(PyCosHarmonic, dh_of_psi_dt);
py_eval_harmonic!(PyCosHarmonic, dh_of_psip_dt);

// ===============================================================================================

#[derive(Clone)]
#[pyclass(name = "_PyNcHarmonic", frozen)]
pub struct PyNcHarmonic(pub(crate) NcHarmonic);

#[pymethods]
impl PyNcHarmonic {
    #[new]
    #[pyo3(signature = (path, interp_type, m, n, phase_method))]
    pub fn new<'py>(
        path: &str,
        interp_type: &str,
        m: i64,
        n: i64,
        phase_method: Option<Bound<'py, PyAny>>,
    ) -> Result<Self, PyEqError> {
        let path = std::path::PathBuf::from(path);
        let method: PhaseMethod = Self::resolve_phase_method(phase_method)?;
        let builder = NcHarmonicBuilder::new(&path, interp_type, m, n).with_phase_method(method);
        Ok(Self(builder.build()?))
    }
}

impl PyNcHarmonic {
    /// Attempts to create a valid [`PhaseMethod`] object from the corresponding [`PyAny`]
    /// Python argument.
    ///
    /// Returns early if a valid string is found that matches one of the simple variants,
    /// otherwise tries to unpack the [`PyAny`] object into a ("Custom", f64) tuple and cast it
    /// into a PhaseMethod::Custom.
    ///
    /// Returns an Error if both string matching and casting fail.
    fn resolve_phase_method<'py>(arg: Option<Bound<'py, PyAny>>) -> PyResult<PhaseMethod> {
        use PhaseMethod::*;

        if arg.is_none() {
            return Ok(PhaseMethod::Average);
        }
        let arg = arg.unwrap();

        match arg.to_string().to_lowercase().as_str() {
            "zero" => return Ok(Zero),
            "average" => return Ok(Average),
            "resonance" => return Ok(Resonance),
            "interpolation" => return Ok(Interpolation),
            _ => (),
        };

        let tuple = arg.cast::<PyTuple>()?;
        let string: String = tuple.get_item(0)?.extract::<String>()?.to_lowercase();
        let value: f64 = tuple.get_item(1)?.extract::<f64>()?;
        match string.as_str() {
            "custom" if value.is_finite() => Ok(PhaseMethod::Custom(value)),
            _ => Err(PyErr::new::<PyTypeError, _>("Invalid 'phase_method'")),
        }
    }
}

py_debug_impl!(PyNcHarmonic);
py_repr_impl!(PyNcHarmonic);
py_get_path!(PyNcHarmonic);
py_get_netcdf_version!(PyNcHarmonic);
py_get_enum_string!(PyNcHarmonic, equilibrium_type);
py_export_getter!(PyNcHarmonic, interp_type, String);
py_export_getter!(PyNcHarmonic, m, i64);
py_export_getter!(PyNcHarmonic, n, i64);
py_get_enum_string!(PyNcHarmonic, phase_method);
py_export_getter!(PyNcHarmonic, phase_average, Option<f64>);
py_export_getter!(PyNcHarmonic, psi_phase_resonance, Option<f64>);
py_export_getter!(PyNcHarmonic, psip_phase_resonance, Option<f64>);
py_export_getter!(PyNcHarmonic, psi_last, Option<f64>);
py_export_getter!(PyNcHarmonic, psip_last, Option<f64>);
py_get_enum_string!(PyNcHarmonic, psi_state);
py_get_enum_string!(PyNcHarmonic, psip_state);
py_get_numpy1D_fallible!(PyNcHarmonic, psi_array);
py_get_numpy1D_fallible!(PyNcHarmonic, psip_array);
py_get_numpy1D!(PyNcHarmonic, alpha_array);
py_get_numpy1D!(PyNcHarmonic, phase_array);
py_eval_harmonic!(PyNcHarmonic, alpha_of_psi);
py_eval_harmonic!(PyNcHarmonic, alpha_of_psip);
py_eval_harmonic!(PyNcHarmonic, phase_of_psi);
py_eval_harmonic!(PyNcHarmonic, phase_of_psip);
py_eval_harmonic!(PyNcHarmonic, h_of_psi);
py_eval_harmonic!(PyNcHarmonic, h_of_psip);
py_eval_harmonic!(PyNcHarmonic, dh_dpsi);
py_eval_harmonic!(PyNcHarmonic, dh_dpsip);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dtheta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dtheta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dzeta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dzeta);
py_eval_harmonic!(PyNcHarmonic, dh_of_psi_dt);
py_eval_harmonic!(PyNcHarmonic, dh_of_psip_dt);

// ===============================================================================================
// ===============================================================================================

/// Unfortunately, since [`Perturbation`] is generic over [`Harmonic`] and pyo3 forbids generics in
/// wrapped types, we must create a new `Py<>Perturbation` object for each [`Harmonic`] implementor.
///
/// Fortunately, all exported methods are identical. However, creating a trait for this is a mess
/// due to pyo3's attributes (if possible at all), so we create this macro instead.
macro_rules! PyPerturbationImpl {
    ($py_perturbation:ident, $py_harmonic:ident, $harmonic:ident) => {
        #[pymethods]
        impl $py_perturbation {
            #[new]
            #[pyo3(signature = (harmonics))]
            pub fn new<'py>(harmonics: Bound<'py, PyList>) -> Result<Self, PyEvalError> {
                let pyharmonics: Vec<$py_harmonic> = harmonics
                    .iter()
                    .map(|h| {
                        h.extract()
                            .expect("Error trying to extract Harmonics from list")
                    })
                    .collect();
                let harmonics: Vec<$harmonic> = pyharmonics.into_iter().map(|h| h.0).collect();
                Ok(Self(Perturbation::new(&harmonics)))
            }

            /// Returns a Python list with all the contained harmonics.
            #[getter]
            pub fn harmonics(&self) -> Vec<$py_harmonic> {
                self.0
                    .harmonics()
                    .iter()
                    .map(|h| $py_harmonic(h.clone()))
                    .collect()
            }

            /// Makes Python type indexable
            pub fn __getitem__(&self, index: usize) -> $py_harmonic {
                $py_harmonic(self.0[index].clone())
            }

            /// Returns the number of the contained harmonics
            pub fn __len__(&self) -> usize {
                self.0.count()
            }
        }

        py_debug_impl!($py_perturbation);
        py_repr_impl!($py_perturbation);
        py_eval_perturbation!($py_perturbation, p_of_psi);
        py_eval_perturbation!($py_perturbation, p_of_psip);
        py_eval_perturbation!($py_perturbation, dp_dpsi);
        py_eval_perturbation!($py_perturbation, dp_dpsip);
        py_eval_perturbation!($py_perturbation, dp_of_psi_dtheta);
        py_eval_perturbation!($py_perturbation, dp_of_psip_dtheta);
        py_eval_perturbation!($py_perturbation, dp_of_psi_dzeta);
        py_eval_perturbation!($py_perturbation, dp_of_psip_dzeta);
        py_eval_perturbation!($py_perturbation, dp_of_psi_dt);
        py_eval_perturbation!($py_perturbation, dp_of_psip_dt);
    };
}

// ===============================================================================================

#[pyclass(name = "_PyCosPerturbation", frozen, dict)]
pub struct PyCosPerturbation(pub(crate) Perturbation<CosHarmonic>);

PyPerturbationImpl!(PyCosPerturbation, PyCosHarmonic, CosHarmonic);

// ===============================================================================================

#[pyclass(name = "_PyNcPerturbation", frozen, dict)]
pub struct PyNcPerturbation(pub(crate) Perturbation<NcHarmonic>);

PyPerturbationImpl!(PyNcPerturbation, PyNcHarmonic, NcHarmonic);
