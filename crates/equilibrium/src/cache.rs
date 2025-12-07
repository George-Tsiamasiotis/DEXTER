//! Helper structs for caching values to avoid unnecessary recalculations.

use rsl_interpolation::Accelerator;

use crate::Harmonic;
use crate::Result;
use crate::harmonics::NcHarmonic;

/// Holds a Harmonic's values evalutated at a specific point.
///
/// Since all the harmonic's methods are called consecutively over the same coordinates, most terms
/// do not need to be calculated every time.
///
/// Similar to the Accelerator, it is an independend object, and does not affect the behavior of the
/// equilibrium objects themselves. It only holds values, calculated by the [`Harmonic`]'s methods.
///
/// ## Note
///
/// The cache should be cloned in each new state calculated from the Stepper.
#[derive(Clone)]
pub struct HarmonicCache {
    hits: usize,
    misses: usize,

    psip: f64,
    theta: f64,
    zeta: f64,

    pub(crate) alpha: f64,
    pub(crate) da_dpsip: f64,
    pub(crate) phase: f64,
    /// The angular part that appears in the sin/cos.
    pub(crate) angle: f64,
    pub(crate) sin: f64,
    pub(crate) cos: f64,
}

impl HarmonicCache {
    /// Creates a new [`HarmonicCache`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Checks if the cache's fields are valid.
    ///
    /// Comparing floats is OK here since they are simply copied between every call, and we
    /// **want** the check to fail with the slightest difference.
    pub(crate) fn is_updated(&mut self, psip: f64, theta: f64, zeta: f64) -> bool {
        if (self.psip == psip) && (self.theta == theta) && (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    /// Updates the cache's fields.
    ///
    /// All the `alpha`, `dalpha`, `phase` values must be calculated as the harmonic itself defines
    /// them.
    pub(crate) fn update(
        &mut self,
        h: &NcHarmonic,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
    ) -> Result<()> {
        self.psip = psip;
        self.theta = theta;
        self.zeta = zeta;
        self.alpha = h.a(psip, acc)?;
        self.phase = h.phase(psip, acc)?;
        self.da_dpsip = h.da_dpsip(psip, acc)?;
        self.angle = h.mod_arg(psip, theta, zeta, acc)?;
        // On some platforms this might be faster than calculating them seperately,
        (self.sin, self.cos) = self.angle.sin_cos();
        Ok(())
    }

    /// Returns the Cache's hit count.
    pub fn hits(&self) -> usize {
        self.hits
    }

    /// Returns the Cache's miss count.
    pub fn misses(&self) -> usize {
        self.misses
    }
}

impl std::fmt::Debug for HarmonicCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HarmonicCache")
            .field("hits  ", &self.hits)
            .field("misses", &self.misses)
            .finish()
    }
}

/// Just in case an initial condition of `ψp=0` actually makes sense.
impl Default for HarmonicCache {
    fn default() -> Self {
        Self {
            hits: Default::default(),
            misses: Default::default(),
            psip: f64::NAN,
            theta: f64::NAN,
            zeta: f64::NAN,
            alpha: f64::NAN,
            phase: f64::NAN,
            da_dpsip: f64::NAN,
            angle: f64::NAN,
            sin: f64::NAN,
            cos: f64::NAN,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::extract::STUB_TEST_NETCDF_PATH;
    use crate::*;
    use std::path::PathBuf;

    #[test]
    fn test_harmonic_cache_update() {
        let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
        let typ = "steffen";
        let m = 2;
        let n = 1;
        let builder = harmonics::NcHarmonicBuilder::new(&path, typ, m, n);
        let harmonic = builder.build().unwrap();
        let mut cache = HarmonicCache::new();
        let mut acc = Accelerator::new();

        assert_eq!(cache.hits(), 0);
        assert_eq!(cache.misses(), 0);

        assert!(cache.psip.is_nan());
        assert!(cache.theta.is_nan());
        assert!(cache.zeta.is_nan());
        assert!(cache.alpha.is_nan());
        assert!(cache.da_dpsip.is_nan());
        assert!(cache.phase.is_nan());
        assert!(cache.angle.is_nan());
        assert!(cache.sin.is_nan());
        assert!(cache.cos.is_nan());

        let psip = 0.002;
        let theta = 3.14;
        let zeta = 3.2;
        assert!(!cache.is_updated(psip, theta, zeta));
        assert_eq!(cache.hits(), 0);
        assert_eq!(cache.misses(), 1);

        assert!(cache.psip.is_nan());
        assert!(cache.theta.is_nan());
        assert!(cache.zeta.is_nan());
        assert!(cache.alpha.is_nan());
        assert!(cache.da_dpsip.is_nan());
        assert!(cache.phase.is_nan());
        assert!(cache.angle.is_nan());
        assert!(cache.sin.is_nan());
        assert!(cache.cos.is_nan());

        cache
            .update(&harmonic, psip, theta, zeta, &mut acc)
            .unwrap();
        assert!(cache.is_updated(psip, theta, zeta));
        assert_eq!(cache.hits(), 1);
        assert_eq!(cache.misses(), 1);

        assert_eq!(cache.psip, psip);
        assert_eq!(cache.theta, theta);
        assert_eq!(cache.zeta, zeta);
        assert!(cache.alpha.is_finite());
        assert!(cache.da_dpsip.is_finite());
        assert!(cache.phase.is_finite());
        assert!(cache.angle.is_finite());
        assert!(cache.sin.is_finite());
        assert!(cache.cos.is_finite());
    }
}
