//! Representation of a total perturbation, a sum of multiple harmonics.

use crate::{CosHarmonic, EvalError, Harmonic};

/// A sum of an arbitrary number of [`Harmonics`](Harmonic).
///
/// It has the general form `Σ{ Φ(n,m)(ψ/ψp, θ, ζ, t)}`, where `Φ(n,m)` the different harmonics.
///
/// Evaluation methods return the sum of the every corresponding evaluation method in each Harmonic.
///
/// The harmonics are stored in the same order as passed on [`Perturbation::new`].
#[non_exhaustive]
pub struct Perturbation<H>(Vec<H>)
where
    H: Harmonic;

impl Perturbation<CosHarmonic> {
    /// Returns a Perturbation without any Harmonics, corresponding to an unperturbed state.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::zero();
    /// # Ok::<_, EqError>(())
    /// ```
    #[must_use]
    pub fn zero() -> Self {
        Self(vec![])
    }
}

impl<H> Perturbation<H>
where
    H: Harmonic,
{
    /// Creates a new Perturbation from an arbitrary number of Harmonics **of the same type**.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// // from analytical harmonics
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    ///
    /// // from numerical harmonics
    /// let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = Perturbation::new(&[
    ///     NcHarmonicBuilder::new(&path, "cubic", 2, 1).build()?,
    ///     NcHarmonicBuilder::new(&path, "cubic", 2, 2).build()?,
    ///     NcHarmonicBuilder::new(&path, "cubic", 3, 1).build()?,
    ///     NcHarmonicBuilder::new(&path, "cubic", 3, 2).build()?,
    /// ]);
    /// # Ok::<_, EqError>(())
    /// ```
    #[must_use]
    pub fn new(harmonics: &[H]) -> Self {
        Self(harmonics.to_vec())
    }

    /// Returns a [`Vec`] of the contained harmonics.
    #[must_use]
    pub fn harmonics(&self) -> Vec<H> {
        self.0.clone()
    }

    /// Returns the number of the contained harmonics.
    #[must_use]
    pub fn count(&self) -> usize {
        self.0.len()
    }

    /// Generates a [`Vec`] of the corresponding caching objects.
    ///
    /// The vec has the same length as the `harmonics` vec, and should be used for evaluations.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches: Vec<CosHarmonicCache> = perturbation.generate_caches();
    /// assert_eq!(perturbation.count(), caches.len());
    ///
    /// let p = perturbation.p_of_psi(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    #[must_use]
    pub fn generate_caches(&self) -> Vec<H::Cache> {
        self.0.iter().map(H::generate_cache).collect()
    }
}

/// Evaluations.
impl<H> Perturbation<H>
where
    H: Harmonic,
{
    /// Calculates the Perturbation's value as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let p_of_psi = perturbation.p_of_psi(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn p_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .h_of_psi(psi, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's value as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let p_of_psip = perturbation.p_of_psip(0.015, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn p_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .h_of_psip(psip, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ψ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_dpsi = perturbation.dp_dpsi(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_dpsi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_dpsi(psi, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ψp`, as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_dpsip = perturbation.dp_dpsip(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_dpsip(psip, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `θ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psi_dtheta = perturbation.dp_of_psi_dtheta(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psi_dtheta(psi, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `θ`, as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psip_dtheta = perturbation.dp_of_psip_dtheta(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psip_dtheta(psip, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ζ`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psi_dzeta = perturbation.dp_of_psi_dzeta(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psi_dzeta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psi_dzeta(psi, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `ζ`, as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psip_dzeta = perturbation.dp_of_psip_dzeta(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psip_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psip_dzeta(psip, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `t`, as a function of `(ψ, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psi_dt = perturbation.dp_of_psi_dt(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psi_dt(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psi_dt(psi, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }

    /// Calculates the Perturbation's derivative with respect to `t`, as a function of `(ψp, θ, ζ, t)`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let perturbation = Perturbation::new(&[
    ///     CosHarmonic::new(1e-3, 1, 2, 0.0),
    ///     CosHarmonic::new(1e-3, 1, 3, 0.0),
    /// ]);
    /// let mut caches = perturbation.generate_caches();
    /// let dp_of_psip_dt = perturbation.dp_of_psip_dt(0.01, 3.14, 3.14, 0.0, &mut caches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the the evaluations fail for any reason.
    pub fn dp_of_psip_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        caches: &mut [H::Cache],
    ) -> Result<f64, EvalError> {
        self.0
            .iter()
            .enumerate()
            .try_fold(0.0, |accumulator, tuple| {
                let (index, harmonic) = tuple;
                harmonic
                    .dh_of_psip_dt(psip, theta, zeta, t, &mut caches[index])
                    .map(|val| accumulator + val)
            })
    }
}

impl<H, Idx> std::ops::Index<Idx> for Perturbation<H>
where
    H: Harmonic,
    Idx: std::slice::SliceIndex<[H], Output = H>,
{
    type Output = H;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.0[index]
    }
}

impl<H> std::fmt::Debug for Perturbation<H>
where
    H: Harmonic + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Perturbation")
            .field("number of harmonics", &self.0.len())
            .finish()
    }
}

#[cfg(test)]
mod perturbation_evals {
    use std::path::PathBuf;

    use crate::extract::TEST_NETCDF_PATH;
    use crate::{HarmonicCache, NcHarmonic, NcHarmonicBuilder};
    use approx::assert_relative_eq;

    use super::*;

    fn create_cos_perturbation() -> Perturbation<CosHarmonic> {
        Perturbation::new(&[
            CosHarmonic::new(1.0, 2, 3, 4.0),
            CosHarmonic::new(5.0, 6, 7, 8.0),
            CosHarmonic::new(9.0, 1, 2, 3.0),
        ])
    }

    #[rustfmt::skip]
    fn create_nc_perturbation() -> Perturbation<NcHarmonic> {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        Perturbation::new(&[
            NcHarmonicBuilder::new(&path, "cubic", 2, 1).build().unwrap(),
            NcHarmonicBuilder::new(&path, "cubic", 2, 1).build().unwrap(),
            NcHarmonicBuilder::new(&path, "cubic", 2, 1).build().unwrap(),
        ])
    }

    #[test]
    fn empty_perturbation() {
        let per = Perturbation::zero();
        let c = &mut per.generate_caches();
        let (p, t) = (0.01, 0.0); // Not used
        assert_eq!(per.p_of_psi(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.p_of_psip(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_dpsi(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_dpsip(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psi_dtheta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psip_dtheta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psi_dzeta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psip_dzeta(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psi_dt(p, 10.0, 20.0, t, c).unwrap(), 0.0);
        assert_eq!(per.dp_of_psip_dt(p, 10.0, 20.0, t, c).unwrap(), 0.0);
    }

    #[test]
    #[rustfmt::skip]
    fn cos_perturbation_evals() {
        let per = create_cos_perturbation();
        let c = &mut per.generate_caches();
        let (p, t) = (0.01, 0.0); // Not used
        let eps = 1e-10;
        assert_relative_eq!(per.p_of_psi(p, 10.0, 20.0, t, c).unwrap(), -7.5934659096, epsilon = eps);
        assert_relative_eq!(per.p_of_psip(p, 10.0, 20.0, t, c).unwrap(), -7.5934659096, epsilon = eps);
        assert_relative_eq!(per.dp_dpsi(p, 10.0, 20.0, t, c).unwrap(), 0.0, epsilon = eps);
        assert_relative_eq!(per.dp_dpsip(p, 10.0, 20.0, t, c).unwrap(), 0.0, epsilon = eps);
        assert_relative_eq!(per.dp_of_psi_dtheta(p, 10.0, 20.0, t, c).unwrap(), 14.2385265316, epsilon = eps);
        assert_relative_eq!(per.dp_of_psip_dtheta(p, 10.0, 20.0, t, c).unwrap(), 14.2385265316, epsilon = eps);
        assert_relative_eq!(per.dp_of_psi_dzeta(p, 10.0, 20.0, t, c).unwrap(), -23.1232478476, epsilon = eps);
        assert_relative_eq!(per.dp_of_psip_dzeta(p, 10.0, 20.0, t, c).unwrap(), -23.1232478476, epsilon = eps);
        assert_relative_eq!(per.dp_of_psi_dt(p, 10.0, 20.0, t, c).unwrap(), 0.0, epsilon = eps);
        assert_relative_eq!(per.dp_of_psip_dt(p, 10.0, 20.0, t, c).unwrap(), 0.0, epsilon = eps);
    }

    #[test]
    #[rustfmt::skip]
    fn nc_perturbation_evals() {
        let per = create_nc_perturbation();
        let c = &mut per.generate_caches();
        let (p, t) = (0.01, 0.0); // Not used
        assert!(per.p_of_psi(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.p_of_psip(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_dpsi(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_dpsip(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psi_dtheta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psip_dtheta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psi_dzeta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psip_dzeta(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psi_dt(p, 10.0, 20.0, t, c).unwrap().is_finite());
        assert!(per.dp_of_psip_dt(p, 10.0, 20.0, t, c).unwrap().is_finite());
    }

    #[test]
    #[allow(unused_results)]
    fn perturbation_cache() {
        let per = create_cos_perturbation();
        let mut c = per.generate_caches();
        let (p, t) = (0.01, 0.0); // Not used

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 0);
            assert_eq!(cache.misses(), 0);
        });

        per.p_of_psi(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 0);
            assert_eq!(cache.misses(), 1);
        });

        per.p_of_psip(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 1);
            assert_eq!(cache.misses(), 1);
        });

        per.dp_of_psi_dtheta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.dp_of_psip_dtheta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.dp_of_psi_dzeta(p, 10.0, 20.0, t, &mut c).unwrap();
        per.dp_of_psip_dzeta(p, 10.0, 20.0, t, &mut c).unwrap();

        c.iter().for_each(|cache| {
            assert_eq!(cache.hits(), 5);
            assert_eq!(cache.misses(), 1);
        });
    }
}
