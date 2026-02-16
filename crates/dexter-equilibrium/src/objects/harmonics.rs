//! Representation of analytical Harmonics.

use crate::eval::HarmonicCache;
use crate::{EquilibriumType, Harmonic, Result};
use crate::{
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl,
};
use std::f64::consts::TAU;

// ===============================================================================================

/// A simple analytical Harmonic of the form `α*cos(mθ-nζ+φ)`, where `α` and `φ` are constants.
///
/// Since this Harmonic is only dependent on `θ` and `ζ`, the rest of the arguments in evaluation
/// methods are ignored.
///
/// Used in pair with [`CosHarmonicCache`].
#[non_exhaustive]
#[derive(Clone)]
pub struct CosHarmonic {
    equilibrium_type: EquilibriumType,
    pub(crate) alpha: f64,
    pub(crate) m: i64,
    pub(crate) n: i64,
    pub(crate) phase: f64,
}

impl CosHarmonic {
    /// Creates a new [`CosHarmonic`].
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let harmonic = CosHarmonic::new(1e-3, 3, 2, 0.0);
    /// ```
    pub fn new(alpha: f64, m: i64, n: i64, phase: f64) -> Self {
        let _m = m as f64;
        let _n = n as f64;
        Self {
            equilibrium_type: EquilibriumType::Analytical,
            alpha,
            m,
            n,
            phase,
        }
    }

    /// Returns the Harmonic's constant amplitude `α`.
    pub fn alpha(&self) -> f64 {
        self.alpha
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
            .field("amplitude", &self.alpha)
            .field("poloidal number `m`", &self.m)
            .field("toroidal number `n`", &self.n)
            .field("phase", &self.phase)
            .finish()
    }
}

/// Stores a [`CosHarmonic`]'s constant parameters and cached quantities
#[derive(Debug)]
pub struct CosHarmonicCache {
    pub(crate) hits: usize,
    pub(crate) misses: usize,
    pub(crate) alpha: f64,
    pub(crate) m: f64,
    pub(crate) n: f64,
    pub(crate) phase: f64,
    /// ====== Order
    /// 0. theta
    /// 1. zeta
    /// 2. modarg
    /// 3. sin
    /// 4. cos
    pub(crate) cache: [f64; 5],
}

impl HarmonicCache for CosHarmonicCache {
    fn is_updated(&mut self, _: f64, theta: f64, zeta: f64, _: f64) -> bool {
        if (self.cache[0] == theta) && (self.cache[1] == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, _: f64, theta: f64, zeta: f64, _: f64) -> Result<()> {
        self.cache[0] = theta;
        self.cache[1] = zeta;
        self.cache[2] = (self.m * theta - self.n * zeta + self.phase).rem_euclid(TAU);
        (self.cache[3], self.cache[4]) = self.cache[2].sin_cos();
        Ok(())
    }
}

harmonic_cache_counts_getter_impl!(CosHarmonicCache);

#[rustfmt::skip]
impl Harmonic for CosHarmonic {
    type Cache = CosHarmonicCache;

    fn generate_cache(&self) -> Self::Cache {
        Self::Cache {
            hits: 0,
            misses: 0,
            alpha: self.alpha,
            m: self.m as f64,
            n: self.n as f64,
            phase: self.phase,
            cache: [f64::NAN; 5],
        }
    }

    fn alpha_of_psi(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.alpha)
    }

    fn alpha_of_psip(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.alpha)
    }

    fn phase_of_psi(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.phase)
    }

    fn phase_of_psip(&self, _: f64, _: f64, _: f64, _: f64, _: &mut CosHarmonicCache) -> Result<f64> {
        Ok(self.phase)
    }

    fn h_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t)?;
        }
        Ok(cache.alpha * cache.cache[4])
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
            cache.update(psi, theta, zeta, t)?;
        }
        Ok(-cache.m * cache.alpha * cache.cache[3])
    }

    fn dh_of_psip_dtheta(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        self.dh_of_psi_dtheta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dzeta(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut CosHarmonicCache) -> Result<f64> {
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t)?;
        }
        Ok(cache.n * cache.alpha * cache.cache[3])
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
        let c = &mut har.generate_cache();

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

        assert!(is_close!(har.alpha_of_psi(p, theta, zeta, t, c)?, 10.0));
        assert!(is_close!(har.alpha_of_psip(p, theta, zeta, t, c)?, 10.0));
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
        let c = &mut h.generate_cache();
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

        h.dh_of_psi_dzeta(psi, theta / 2.0, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 2);
        assert_eq!(c.misses(), 3);

        // dbg!(&c); // Accelerator counts should also match
    }
}
