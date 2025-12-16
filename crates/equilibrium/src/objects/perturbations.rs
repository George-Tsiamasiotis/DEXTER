//! Representation of a total perturbation, a sum of multiple harmonics.

use rsl_interpolation::Accelerator;

use crate::NcHarmonic;
use crate::cache::HarmonicCache;
use crate::{Harmonic, Perturbation};

/// A sum of different perturbation [`NcHarmonics`](NcHarmonic).
///
/// It has the general form
///     `Σ{ α(n,m)(ψp) * cos(mθ-nζ+φ0) }`.
pub struct NcPerturbation {
    harmonics: Vec<NcHarmonic>,
}

// Creation and data extraction
impl NcPerturbation {
    /// Creates a Perturbation from different [`NcHarmonics`](NcHarmonic).
    ///
    /// # Examples
    ///
    /// No perturbations:
    /// ```
    /// # use equilibrium::*;
    /// let perturbation = NcPerturbation::from_harmonics(&[]);
    /// ```
    ///
    /// Multiple perturbations:
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # let path = PathBuf::from(extract::STUB_TEST_NETCDF_PATH);
    /// let perturbation = NcPerturbation::from_harmonics(&[
    ///     NcHarmonicBuilder::new(&path, "steffen", 2, 1).build()?,
    ///     NcHarmonicBuilder::new(&path, "steffen", 3, 2).build()?,
    /// ]);
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn from_harmonics(harmonics: &[NcHarmonic]) -> Self {
        Self {
            harmonics: harmonics.into(),
        }
    }

    pub fn get_harmonics(&self) -> Vec<NcHarmonic> {
        self.harmonics.clone()
    }
}

impl Perturbation for NcPerturbation {
    fn p(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .h(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dpsip(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dtheta(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dzeta(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn dp_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> crate::Result<f64> {
        self.harmonics.iter().enumerate().try_fold(0.0, |p, tuple| {
            let (index, harmonic) = tuple;
            harmonic
                .dh_dt(psip, theta, zeta, acc, &mut caches[index])
                .map(|v| p + v)
        })
    }

    fn len(&self) -> usize {
        self.harmonics.len()
    }
}

impl std::fmt::Debug for NcPerturbation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for harmonic in self.get_harmonics() {
            let _ = harmonic.fmt(f);
        }
        Ok(())
    }
}
