//! Representation of A single analytical [`Harmonic`] of the form `Vcos(mθ-nζ+φ)`.

use crate::eval::HarmonicCache;
use crate::{EquilibriumType, Harmonic, Result};
use crate::{
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl,
};
use std::f64::consts::TAU;

/// A simple analytical Harmonic of the form `α*cos(mθ-nζ+φ)`, where `α` and `φ` are constants.
///
/// Since this Harmonic is only dependent on `θ` and `ζ`, the rest of the arguments in evaluation
/// methods are ignored.
///
/// Used in pair with [`CosHarmonicCache`].
#[non_exhaustive]
pub struct CosHarmonic {
    equilibrium_type: EquilibriumType,
    pub(crate) ampl: f64,
    pub(crate) m: i64,
    pub(crate) n: i64,
    pub(crate) phase: f64,

    // To be used in actual calculations.
    _m: f64,
    _n: f64,

    // [m, n, phase]
    _cache_args: [f64; 3],
}

impl CosHarmonic {
    /// Creates a new [`CosHarmonic`].
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let harmonic = CosHarmonic::new(1e-3, 3, 2, 0.0);
    /// ```
    pub fn new(ampl: f64, m: i64, n: i64, phase: f64) -> Self {
        let _m = m as f64;
        let _n = n as f64;
        Self {
            equilibrium_type: EquilibriumType::Analytical,
            ampl,
            m,
            n,
            phase,
            _m,
            _n,
            _cache_args: [_m, _n, phase],
        }
    }

    /// Returns the Harmonic's constant amplitude `α`.
    pub fn ampl(&self) -> f64 {
        self.ampl
    }

    /// Returns the Harmonic's constant phase `φ`.
    pub fn phase(&self) -> f64 {
        self.phase
    }

    harmonic_mode_number_getter_impl!();
    equilibrium_type_getter_impl!();
}

impl std::fmt::Debug for CosHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CosHarmonic")
            .field("equilibrium_type", &self.equilibrium_type)
            .field("amplitude", &self.ampl)
            .field("poloidal number `m`", &self.m)
            .field("toroidal number `n`", &self.n)
            .field("phase", &self.phase)
            .finish()
    }
}

/// Cache for [`CosHarmonic`] values.
///
/// # Note
///
/// Do not use the same NcHarmonicCache for both `ψ` and `ψp` evaluations, since this would
/// invalidate the cached values.
///
/// # Example
///
/// ```
/// # use std::path::PathBuf;
/// # use dexter_equilibrium::*;
/// let harmonic = CosHarmonic::new(1e-3, 3, 2, 0.0);
/// let mut cache = harmonic.get_default_cache();
/// let h = harmonic.h_of_psi(0.0, 0.2, 0.3, 0.0, &mut cache)?;
/// # Ok::<_, EqError>(())
/// ```
///
#[derive(Debug)]
#[non_exhaustive]
pub struct CosHarmonicCache {
    hits: usize,
    misses: usize,
    theta: f64,
    zeta: f64,
    modarg: f64,
    cos: f64,
    sin: f64,
}

harmonic_cache_counts_getter_impl!(CosHarmonicCache);

impl HarmonicCache for CosHarmonicCache {
    fn is_updated(&mut self, _: f64, theta: f64, zeta: f64, _: f64) -> bool {
        if (self.theta == theta) && (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, _: f64, theta: f64, zeta: f64, _: f64, args: &[f64]) -> Result<()> {
        let m = args.first().expect("Always exists");
        let n = args.get(1).expect("Always exists");
        let phase = args.get(2).expect("Always exists");

        self.theta = theta;
        self.zeta = zeta;
        self.modarg = (m * theta - n * zeta + phase).rem_euclid(TAU);
        (self.sin, self.cos) = self.modarg.sin_cos();
        Ok(())
    }

    fn sin(&self) -> f64 {
        self.sin
    }

    fn cos(&self) -> f64 {
        self.cos
    }

    fn ampl(&self) -> f64 {
        unreachable!("return the constant amplitude")
    }

    fn set_ampl(&mut self, _: f64) {
        unreachable!("ampl is constant")
    }

    fn set_phase(&mut self, _: f64) {
        unreachable!("phase is constant")
    }

    fn set_dampl(&mut self, _: f64) {
        unreachable!("ampl is constant")
    }

    fn phase(&self) -> f64 {
        unreachable!("phase not stored in cache")
    }

    fn dampl(&self) -> f64 {
        unreachable!("ampl is constant")
    }

    fn flux_acc(&mut self) -> &mut rsl_interpolation::Accelerator {
        unreachable!("Analytical object has no Accelerator")
    }
}

impl Default for CosHarmonicCache {
    /// Set all initial values to NaN so first is_updated always returns false.
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            theta: f64::NAN,
            zeta: f64::NAN,
            modarg: f64::NAN,
            cos: f64::NAN,
            sin: f64::NAN,
        }
    }
}

#[rustfmt::skip]
impl Harmonic for CosHarmonic {
    type Cache = CosHarmonicCache;

    fn ampl_of_psi(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.ampl)
    }

    fn ampl_of_psip(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.ampl)
    }

    fn phase_of_psi(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.phase)
    }

    fn phase_of_psip(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.phase)
    }

    fn h_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t, &self._cache_args)?;
        }
        Ok(self.ampl * cache.cos())
    }

    fn h_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        self.h_of_psi(psip, theta, zeta, t, cache)
    }

    fn dh_dpsi(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }

    fn dh_dpsip(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }

    fn dh_of_psi_dtheta(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t, &self._cache_args)?;
        }
        Ok(-self._m * self.ampl * cache.sin())
    }

    fn dh_of_psip_dtheta(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        self.dh_of_psi_dtheta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dzeta(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t, &self._cache_args)?;
        }
        Ok(self._n * self.ampl * cache.sin())
    }

    fn dh_of_psip_dzeta(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        self.dh_of_psi_dzeta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dt(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }

    fn dh_of_psip_dt(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }
}

#[cfg(test)]
mod cos_harmonic_values {
    use super::*;
    use is_close::is_close;

    #[test]
    #[rustfmt::skip]
    fn desmos_values() -> Result<()> {
        let har = dbg!(CosHarmonic::new(10.0, 3, 2, 1.0));
        let c = &mut har.get_default_cache();

        let p = 0.0; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        assert!(is_close!(har.h_of_psi(p, theta, zeta, t, c)?, 5.40302305868));
        assert!(is_close!(har.h_of_psip(p, theta, zeta, t, c)?, 5.40302305868));
        assert!(is_close!(har.dh_of_psi_dtheta(p, theta, zeta, t, c)?, -25.2441295442));
        assert!(is_close!(har.dh_of_psip_dtheta(p, theta, zeta, t, c)?, -25.2441295442));
        assert!(is_close!(har.dh_of_psi_dzeta(p, theta, zeta, t, c)?, 16.8294196962));
        assert!(is_close!(har.dh_of_psip_dzeta(p, theta, zeta, t, c)?, 16.8294196962));

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 5);

        assert!(is_close!(har.ampl_of_psi(p, theta, zeta, t, c)?, 10.0));
        assert!(is_close!(har.ampl_of_psip(p, theta, zeta, t, c)?, 10.0));
        assert!(is_close!(har.phase_of_psi(p, theta, zeta, t, c)?, 1.0));
        assert!(is_close!(har.phase_of_psip(p, theta, zeta, t, c)?, 1.0));
        assert!(is_close!(har.dh_dpsi(p, theta, zeta, t, c)?, 0.0));
        assert!(is_close!(har.dh_dpsip(p, theta, zeta, t, c)?, 0.0));
        assert!(is_close!(har.dh_of_psi_dt(p, theta, zeta, t, c)?, 0.0));
        assert!(is_close!(har.dh_of_psip_dt(p, theta, zeta, t, c)?, 0.0));
        Ok(())
    }
}

#[cfg(test)]
mod cos_harmonic_cache {

    use super::*;

    #[test]
    fn counts() {
        let h = dbg!(CosHarmonic::new(10.0, 3, 2, 1.0));
        let c = &mut h.get_default_cache();
        let psi = 0.01; // not checked
        let t = 0.0; // not checked

        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        h.phase_of_psi(psi, 0.1, 0.1, t, c).unwrap(); // Does not check
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        h.dh_of_psi_dtheta(0.01, 0.1, 0.1, t, c).unwrap(); // First check
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 1);

        let theta = 3.14;
        let zeta = 1.0;
        h.h_of_psi(psi, theta, zeta, t, c).unwrap();
        h.dh_of_psi_dtheta(psi, theta, zeta, t, c).unwrap();
        h.dh_of_psi_dzeta(psi, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 2);
        assert_eq!(c.misses(), 2);
        // dbg!(&c); // Accelerator counts should also match
    }
}
