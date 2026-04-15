//! `dexter_equilibrium::Perturbation` newtype, constructor and method exports.

use dexter::dexter_equilibrium::{CosHarmonic, NcHarmonic, Perturbation};
use pyo3::prelude::*;
use pyo3::types::PyList;

use crate::pyerror::PyEvalError;
use crate::pyharmonics::{PyCosHarmonic, PyNcHarmonic};
use crate::{py_debug_impl, py_eval_perturbation, py_repr_impl};

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
