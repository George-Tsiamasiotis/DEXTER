//! `dexter_simulate::COMs` newtype, constructors and method exports.

use dexter::dexter_simulate::COMs;
use ndarray::Array1;
use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::{py_debug_impl, py_export_pub_field, py_repr_impl};

// ===============================================================================================

#[pyclass(name = "_PyCOMs", frozen)]
pub struct PyCOMs(pub(crate) COMs);

#[pymethods]
impl PyCOMs {
    #[new]
    #[pyo3(signature = (/, energy, pzeta, mu))]
    pub fn new(energy: Option<f64>, pzeta: Option<f64>, mu: Option<f64>) -> Self {
        Self(COMs { energy, pzeta, mu })
    }
}

py_debug_impl!(PyCOMs);
py_repr_impl!(PyCOMs);
py_export_pub_field!(PyCOMs, energy, Option<f64>);
py_export_pub_field!(PyCOMs, pzeta, Option<f64>);
py_export_pub_field!(PyCOMs, mu, Option<f64>);

// ===============================================================================================

/// See [`py_particle_generics_impl`].
#[allow(non_camel_case_types)]
mod py_coms_generics_impl {
    use super::*;

    use crate::{generic_energy_of_psi_grid_impl, generic_energy_of_psip_grid_impl};

    type uniQ = crate::pyqfactors::PyUnityQfactor;
    type parQ = crate::pyqfactors::PyParabolicQfactor;
    type ncdQ = crate::pyqfactors::PyNcQfactor;

    type larC = crate::pycurrents::PyLarCurrent;
    type ncdC = crate::pycurrents::PyNcCurrent;

    type larB = crate::pybfields::PyLarBfield;
    type ncdB = crate::pybfields::PyNcBfield;

    // ================= `energy_of_psi_grid`

    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_uniQ_larC_larB, uniQ, larC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_uniQ_larC_ncdB, uniQ, larC, ncdB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_uniQ_ncdC_larB, uniQ, ncdC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_uniQ_ncdC_ncdB, uniQ, ncdC, ncdB);

    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_parQ_larC_larB, parQ, larC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_parQ_larC_ncdB, parQ, larC, ncdB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_parQ_ncdC_larB, parQ, ncdC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_parQ_ncdC_ncdB, parQ, ncdC, ncdB);

    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_ncdQ_larC_larB, ncdQ, larC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_ncdQ_larC_ncdB, ncdQ, larC, ncdB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_ncdQ_ncdC_larB, ncdQ, ncdC, larB);
    generic_energy_of_psi_grid_impl!(__energy_of_psi_grid_ncdQ_ncdC_ncdB, ncdQ, ncdC, ncdB);

    // ================= `energy_of_psip_grid`

    generic_energy_of_psip_grid_impl!(__energy_of_psip_grid_larC_larB, larC, larB);
    generic_energy_of_psip_grid_impl!(__energy_of_psip_grid_larC_ncdB, larC, ncdB);
    generic_energy_of_psip_grid_impl!(__energy_of_psip_grid_ncdC_larB, ncdC, larB);
    generic_energy_of_psip_grid_impl!(__energy_of_psip_grid_ncdC_ncdB, ncdC, ncdB);
}
