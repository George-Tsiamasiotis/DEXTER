//! Representation of analytical Harmonics.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl,
};
use std::f64::consts::TAU;

use crate::EvalError;
use crate::{EquilibriumType, Harmonic, HarmonicCache, LastClosedFluxSurface};

// ===============================================================================================

/// A simple analytical Harmonic of the form `ε*sqrt(ψ/ψlast)*cos(mθ-nζ+φ)` or
/// `ε*sqrt(ψp/ψplast)*cos(mθ-nζ+φ)`, where `ε` and `φ` are constants and `ψlast/ψplast` the last
/// closed flux surface.
///
/// The square root is necessary since the harmonic must behave like the square root of the
/// magnetic flux close to the axis.
///
/// Used in pair with [`CosHarmonicCache`].
#[non_exhaustive]
#[derive(Clone)]
pub struct CosHarmonic {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The harmonic's "amplitude" `ε`. Corresponds the value of the amplitude at the last closed
    /// flux surface.
    epsilon: f64,
    /// The harmonic's poloidal mode number `m`.
    m: i64,
    /// The harmonic's toroidal mode number `n`.
    n: i64,
    /// The harmonic's phase.
    phase: f64,
    /// The last closed flux surface.
    lcfs: LastClosedFluxSurface,
    /// The value of the last closed toroidal flux surface, if the harmonic was defined through the
    /// toroidal flux.
    psi_last: Option<f64>,
    /// The value of the last closed poloidal flux surface, if the harmonic was defined through the
    /// poloidal flux.
    psip_last: Option<f64>,
}

impl CosHarmonic {
    /// Creates a new [`CosHarmonic`].
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    /// let harmonic = CosHarmonic::new(1e-3, lcfs, 3, 2, 0.0);
    /// ```
    #[must_use]
    pub fn new(epsilon: f64, lcfs: LastClosedFluxSurface, m: i64, n: i64, phase: f64) -> Self {
        let psi_last: Option<f64>;
        let psip_last: Option<f64>;
        match lcfs {
            LastClosedFluxSurface::Toroidal(last) => {
                psi_last = Some(last);
                psip_last = None;
            }
            LastClosedFluxSurface::Poloidal(last) => {
                psi_last = None;
                psip_last = Some(last)
            }
        }

        Self {
            equilibrium_type: EquilibriumType::Analytical,
            lcfs,
            epsilon,
            m,
            n,
            phase,
            psi_last,
            psip_last,
        }
    }

    /// Returns the Harmonic's last closed flux surface.
    #[must_use]
    pub fn lcfs(&self) -> LastClosedFluxSurface {
        self.lcfs
    }

    /// Returns the Harmonic's constant "amplitude" `ε`.
    #[must_use]
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    /// Returns the Harmonic's constant phase `φ`.
    #[must_use]
    pub fn phase(&self) -> f64 {
        self.phase
    }

    /// Returns the value of the last closed toroidal flux `ψ_last`.
    #[must_use]
    pub fn psi_last(&self) -> Option<f64> {
        self.psi_last
    }

    /// Returns the value of the last closed poloidal flux `ψp_last`.
    #[must_use]
    pub fn psip_last(&self) -> Option<f64> {
        self.psip_last
    }

    harmonic_mode_number_getter_impl!();
    equilibrium_type_getter_impl!();
}

impl std::fmt::Debug for CosHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CosHarmonic")
            .field("equilibrium_type", &self.equilibrium_type)
            .field("epsilon", &self.epsilon)
            .field("LCFS", &self.lcfs())
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
    /// The harmonic's constant `ε`.
    epsilon: f64,
    /// The harmonic's poloidal mode number `m`, casted to `f64` to use for evaluations.
    m: f64,
    /// The harmonic's toroidal mode number `n`, casted to `f64` to use for evaluations.
    n: f64,
    /// The harmonic's phase.
    phase: f64,
    /// The square root of the value of the last closed flux surface.
    lcfs_root: f64,
    /// Ordered array with the intermediate values.
    ///
    /// cache = [
    ///     0 = flux,
    ///     1 = theta,
    ///     2 = zeta,
    ///     3 = sqrt(flux)
    ///     4 = modarg
    ///     5 = sin,
    ///     6 = cos
    /// ].
    cache: [f64; 7],
}

impl Default for CosHarmonicCache {
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            epsilon: f64::NAN,
            m: f64::NAN,
            n: f64::NAN,
            phase: f64::NAN,
            lcfs_root: f64::NAN,
            cache: [f64::NAN; 7],
        }
    }
}

impl HarmonicCache for CosHarmonicCache {
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) -> bool {
        #[expect(
            clippy::float_cmp,
            reason = "we want a cache hit only if all values are exactly equal"
        )]
        if (self.cache[0] == flux) && (self.cache[1] == theta) && (self.cache[2] == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) {
        self.cache[0] = flux;
        self.cache[1] = theta;
        self.cache[2] = zeta;
        self.cache[3] = flux.sqrt();
        self.cache[4] = (self.m * theta - self.n * zeta + self.phase).rem_euclid(TAU);
        (self.cache[5], self.cache[6]) = self.cache[4].sin_cos();
    }

    harmonic_cache_counts_getter_impl!(CosHarmonicCache);
}

impl Harmonic for CosHarmonic {
    type Cache = CosHarmonicCache;

    fn generate_cache(&self) -> Self::Cache {
        let lcfs_root = match self.lcfs {
            LastClosedFluxSurface::Toroidal(last) | LastClosedFluxSurface::Poloidal(last) => {
                last.sqrt()
            }
        };
        Self::Cache {
            epsilon: self.epsilon,
            m: self.m as f64,
            n: self.n as f64,
            phase: self.phase,
            lcfs_root,
            ..Default::default()
        }
    }

    fn alpha_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("α(ψ)".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                self.epsilon * cache.cache[3] / cache.lcfs_root
            ))
        }
    }

    fn alpha_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("α(ψp)".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                self.epsilon * cache.cache[3] / cache.lcfs_root
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("h(ψ)".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[6]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("h(ψp)".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[6]
            ))
        }
    }

    fn dh_dpsi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dψ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.epsilon / (2.0 * cache.lcfs_root * cache.cache[3]) * cache.cache[6]
            ))
        }
    }

    fn dh_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut CosHarmonicCache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dψp".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.epsilon / (2.0 * cache.lcfs_root * cache.cache[3]) * cache.cache[6]
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dθ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                -cache.m * cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[5]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dθ".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                -cache.m * cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[5]
            ))
        }
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
        if self.psi_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψ)/dζ".into()))
        } else {
            if !cache.is_updated(psi, theta, zeta, t) {
                cache.update(psi, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.n * cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[5]
            ))
        }
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
        if self.psip_last.is_none() {
            Err(EvalError::UndefinedEvaluation("dh(ψp)/dζ".into()))
        } else {
            if !cache.is_updated(psip, theta, zeta, t) {
                cache.update(psip, theta, zeta, t);
            }
            Ok(debug_assert_is_finite!(
                cache.n * cache.epsilon * cache.cache[3] / cache.lcfs_root * cache.cache[5]
            ))
        }
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
    fn desmos_values_toroidal_lcfs() -> Result<(), EvalError> {
        let lcfs = LastClosedFluxSurface::Toroidal(0.45);
        let har = dbg!(CosHarmonic::new(10.0, lcfs, 3, 2, 1.0));
        let c = &mut har.generate_cache();

        let p = 0.1; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        let eps = 1e-10;
        assert_relative_eq!(har.alpha_of_psi(p, theta, zeta, t, c)?, 4.71404520791, epsilon = eps);
        assert_relative_eq!(har.phase_of_psi(p, theta, zeta, t, c)?, 1.0, epsilon = eps);
        assert_relative_eq!(har.h_of_psi(p, theta, zeta, t, c)?, 2.5470094958, epsilon = eps);
        assert_relative_eq!(har.dh_dpsi(p, theta, zeta, t, c)?, 12.735047479, epsilon = eps);
        assert_relative_eq!(har.dh_of_psi_dtheta(p, theta, zeta, t, c)?, -11.9001967906, epsilon = eps);
        assert_relative_eq!(har.dh_of_psi_dzeta(p, theta, zeta, t, c)?, 7.93346452706, epsilon = eps);

        assert!(har.alpha_of_psip(p, theta, zeta, t,c).is_err());
        assert!(har.h_of_psip(p, theta, zeta, t,c).is_err());
        assert!(har.dh_dpsip(p, theta, zeta, t,c).is_err());
        assert!(har.dh_of_psip_dtheta(p, theta, zeta, t,c).is_err());
        assert!(har.dh_of_psip_dzeta(p, theta, zeta, t,c).is_err());

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 4);

        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn desmos_values_poloidal_lcfs() -> Result<(), EvalError> {
        let lcfs = LastClosedFluxSurface::Poloidal(0.45);
        let har = dbg!(CosHarmonic::new(10.0, lcfs, 3, 2, 1.0));
        let c = &mut har.generate_cache();

        let p = 0.1; // not used
        let theta = 0.2;
        let zeta = 0.3;
        let t = 0.0; // not used

        let eps = 1e-10;
        assert_relative_eq!(har.alpha_of_psip(p, theta, zeta, t, c)?, 4.71404520791, epsilon = eps);
        assert_relative_eq!(har.phase_of_psip(p, theta, zeta, t, c)?, 1.0, epsilon = eps);
        assert_relative_eq!(har.h_of_psip(p, theta, zeta, t, c)?, 2.5470094958, epsilon = eps);
        assert_relative_eq!(har.dh_dpsip(p, theta, zeta, t, c)?, 12.735047479, epsilon = eps);
        assert_relative_eq!(har.dh_of_psip_dtheta(p, theta, zeta, t, c)?, -11.9001967906, epsilon = eps);
        assert_relative_eq!(har.dh_of_psip_dzeta(p, theta, zeta, t, c)?, 7.93346452706, epsilon = eps);

        assert!(har.alpha_of_psi(p, theta, zeta, t,c).is_err());
        assert!(har.h_of_psi(p, theta, zeta, t,c).is_err());
        assert!(har.dh_dpsi(p, theta, zeta, t,c).is_err());
        assert!(har.dh_of_psi_dtheta(p, theta, zeta, t,c).is_err());
        assert!(har.dh_of_psi_dzeta(p, theta, zeta, t,c).is_err());

        assert_eq!(c.misses(), 1);
        assert_eq!(c.hits(), 4);

        Ok(())
    }
}

#[cfg(test)]
#[expect(unused_results)]
mod cos_harmonic_cache {

    use super::*;

    #[test]
    fn counts() {
        let lcfs = LastClosedFluxSurface::Toroidal(0.45);
        let h = dbg!(CosHarmonic::new(10.0, lcfs, 3, 2, 1.0));
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
