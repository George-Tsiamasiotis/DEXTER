//! Representation of analytical Harmonics.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl,
};
use std::f64::consts::TAU;

use crate::EvalError;
use crate::{EquilibriumType, Harmonic, HarmonicCache};

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
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The harmonic's amplitude.
    alpha: f64,
    /// The harmonic's poloidal mode number `m`.
    m: i64,
    /// The harmonic's toroidal mode number `n`.
    n: i64,
    /// The harmonic's phase.
    phase: f64,
}

impl CosHarmonic {
    /// Creates a new [`CosHarmonic`].
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let harmonic = CosHarmonic::new(1e-3, 3, 2, 0.0);
    /// ```
    #[must_use]
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
    #[must_use]
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Returns the Harmonic's constant phase `φ`.
    #[must_use]
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

/// Stores a [`CosHarmonic`]'s constant parameters and cached quantities.
#[derive(Debug, Clone)]
pub struct CosHarmonicCache {
    /// The number of cache hits.
    hits: usize,
    /// The number of cache misses.
    misses: usize,
    /// The harmonic's amplitude.
    alpha: f64,
    /// The harmonic's poloidal mode number `m`, casted to `f64` to use for evaluations.
    m: f64,
    /// The harmonic's toroidal mode number `n`, casted to `f64` to use for evaluations.
    n: f64,
    /// The harmonic's phase.
    phase: f64,
    /// Ordered array with the intermediate values.
    ///
    /// cache = [
    ///     0 = theta,
    ///     1 = zeta,
    ///     2 = modarg,
    ///     3 = sin,
    ///     4 = cos
    /// ].
    cache: [f64; 5],
}

impl Default for CosHarmonicCache {
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            alpha: f64::NAN,
            m: f64::NAN,
            n: f64::NAN,
            phase: f64::NAN,
            cache: [f64::NAN; 5],
        }
    }
}

impl HarmonicCache for CosHarmonicCache {
    fn is_updated(&mut self, _: f64, theta: f64, zeta: f64, _: f64) -> bool {
        #[expect(
            clippy::float_cmp,
            reason = "we want a cache hit only if all values are exactly equal"
        )]
        if (self.cache[0] == theta) && (self.cache[1] == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, _: f64, theta: f64, zeta: f64, _: f64) {
        self.cache[0] = theta;
        self.cache[1] = zeta;
        self.cache[2] = (self.m * theta - self.n * zeta + self.phase).rem_euclid(TAU);
        (self.cache[3], self.cache[4]) = self.cache[2].sin_cos();
    }

    harmonic_cache_counts_getter_impl!(CosHarmonicCache);
}

impl Harmonic for CosHarmonic {
    type Cache = CosHarmonicCache;

    fn generate_cache(&self) -> Self::Cache {
        Self::Cache {
            alpha: self.alpha,
            m: self.m as f64,
            n: self.n as f64,
            phase: self.phase,
            ..Default::default()
        }
    }

    fn alpha_of_psi(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(self.alpha))
    }

    fn alpha_of_psip(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(debug_assert_is_finite!(self.alpha))
    }

    fn phase_of_psi(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(self.phase))
    }

    fn phase_of_psip(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(debug_assert_is_finite!(self.phase))
    }

    fn h_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t);
        }
        Ok(debug_assert_is_finite!(cache.alpha * cache.cache[4]))
    }

    fn h_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        self.h_of_psi(psip, theta, zeta, t, cache)
    }

    fn dh_dpsi(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(0.0_f64))
    }

    fn dh_dpsip(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(debug_assert_is_finite!(0.0_f64))
    }

    fn dh_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t);
        }
        Ok(debug_assert_is_finite!(
            -cache.m * cache.alpha * cache.cache[3]
        ))
    }

    fn dh_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        self.dh_of_psi_dtheta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dzeta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if !cache.is_updated(psi, theta, zeta, t) {
            cache.update(psi, theta, zeta, t);
        }
        Ok(debug_assert_is_finite!(
            cache.n * cache.alpha * cache.cache[3]
        ))
    }

    fn dh_of_psip_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        self.dh_of_psi_dzeta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dt(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(0.0_f64))
    }

    fn dh_of_psip_dt(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(debug_assert_is_finite!(0.0_f64))
    }
}

#[cfg(test)]
mod cos_harmonic_values {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    #[rustfmt::skip]
    fn desmos_values() -> Result<(), EvalError> {
        let har = dbg!(CosHarmonic::new(10.0, 3, 2, 1.0));
        let c = &mut har.generate_cache();

        let p = 0.0; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        let eps = 1e-10;
        assert_relative_eq!(har.h_of_psi(p, theta, zeta, t, c)?, 5.40302305868, epsilon = eps);
        assert_relative_eq!(har.h_of_psip(p, theta, zeta, t, c)?, 5.40302305868, epsilon = eps);
        assert_relative_eq!(har.dh_of_psi_dtheta(p, theta, zeta, t, c)?, -25.2441295442, epsilon = eps);
        assert_relative_eq!(har.dh_of_psip_dtheta(p, theta, zeta, t, c)?, -25.2441295442, epsilon = eps);
        assert_relative_eq!(har.dh_of_psi_dzeta(p, theta, zeta, t, c)?, 16.8294196962, epsilon = eps);
        assert_relative_eq!(har.dh_of_psip_dzeta(p, theta, zeta, t, c)?, 16.8294196962, epsilon = eps);

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 5);

        assert_eq!(har.alpha_of_psi(p, theta, zeta, t, c)?, 10.0);
        assert_eq!(har.alpha_of_psip(p, theta, zeta, t, c)?, 10.0);
        assert_eq!(har.phase_of_psi(p, theta, zeta, t, c)?, 1.0);
        assert_eq!(har.phase_of_psip(p, theta, zeta, t, c)?, 1.0);
        assert_eq!(har.dh_dpsi(p, theta, zeta, t, c)?, 0.0);
        assert_eq!(har.dh_dpsip(p, theta, zeta, t, c)?, 0.0);
        assert_eq!(har.dh_of_psi_dt(p, theta, zeta, t, c)?, 0.0);
        assert_eq!(har.dh_of_psip_dt(p, theta, zeta, t, c)?, 0.0);
        Ok(())
    }
}

#[cfg(test)]
#[expect(unused_results)]
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
    }
}
