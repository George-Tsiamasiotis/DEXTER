//! Helper structs for caching values to avoid unnecessary recalculations.

use std::f64::consts::TAU;

use rsl_interpolation::Accelerator;

use crate::Harmonic;
use crate::Result;
use crate::{Flux, Radians};
use crate::{NcHarmonic, PhaseMethod};

/// Required methods for caching an `impl Harmonic` object's values.
pub trait HarmonicCache {
    /// Returns the Cache's hit count.
    fn hits(&self) -> usize;

    /// Returns the Cache's miss count.
    fn misses(&self) -> usize;

    /// Checks if the cache's fields are valid.
    ///
    /// Comparing floats is OK here since they are simply copied between every call, and we
    /// **want** the check to fail with the slightest difference.
    fn is_updated(&mut self, psip: Flux, theta: Radians, zeta: Radians) -> bool;

    /// Updates the cache's fields.
    fn update(
        &mut self,
        h: &NcHarmonic,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        acc: &mut Accelerator,
    ) -> Result<()>;

    fn alpha(&self) -> f64;

    fn dalpha(&self) -> f64;

    fn cos(&self) -> f64;

    fn sin(&self) -> f64;
}

/// Holds an [`NcHarmonic`]'s values evalutated at a specific point.
///
/// Since all the harmonic's methods are called consecutively over the same coordinates, most terms
/// do not need to be calculated every time.
///
/// Similar to the Accelerators, they are stored inside State, and do not affect the behavior of the
/// equilibrium objects themselves.
///
/// The cache should be cloned in each new state calculated from the Stepper.
#[derive(Clone, Default)]
pub struct NcHarmonicCache {
    hits: usize,
    misses: usize,
    psip: Flux,
    theta: Radians,
    zeta: Radians,
    pub(crate) alpha: f64,
    pub(crate) phase: Radians,
    pub(crate) dalpha: f64,
    pub(crate) sin: f64,
    pub(crate) cos: f64,
}

impl NcHarmonicCache {
    /// Creates a new [`NcHarmonicCache`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the phase value `φ(ψp)`, depending on the harmonic's [`PhaseMethod`].
    fn calculate_phase(h: &NcHarmonic, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        // Options are always Some when the correct method is set
        match h.phase_method {
            PhaseMethod::Zero => Ok(0.0),
            PhaseMethod::Average => Ok(h.phase_average.expect("is Some")),
            PhaseMethod::Resonance => Ok(h.phase_average.expect("is Some")),
            PhaseMethod::Interpolation => Ok(h.phase(psip, acc)?),
            PhaseMethod::Custom(value) => Ok(value),
        }
    }
}

impl HarmonicCache for NcHarmonicCache {
    fn hits(&self) -> usize {
        self.hits
    }

    fn misses(&self) -> usize {
        self.misses
    }

    /// Checks if the cache's fields are valid.
    ///
    /// Comparing floats is OK here since they are simply copied between every call, and we
    /// **want** the check to fail with the slightest difference.
    fn is_updated(&mut self, psip: Flux, theta: Radians, zeta: Radians) -> bool {
        if (self.psip == psip) && (self.theta == theta) && (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    /// Updates the cache's fields.
    fn update(
        &mut self,
        h: &NcHarmonic,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        acc: &mut Accelerator,
    ) -> Result<()> {
        self.psip = psip;
        self.theta = theta;
        self.zeta = zeta;
        self.alpha = h.a(psip, acc)?;
        self.phase = Self::calculate_phase(h, psip, acc)?;
        self.dalpha = h.da_dpsip(psip, acc)?;
        let mod_arg = (h._m * self.theta - h._n * self.zeta + self.phase).rem_euclid(TAU);
        (self.sin, self.cos) = mod_arg.sin_cos();
        Ok(())
    }

    fn alpha(&self) -> f64 {
        self.alpha
    }

    fn dalpha(&self) -> f64 {
        self.dalpha
    }

    fn cos(&self) -> f64 {
        self.cos
    }

    fn sin(&self) -> f64 {
        self.sin
    }
}

impl std::fmt::Debug for NcHarmonicCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HarmonicCache")
            .field("hits  ", &self.hits)
            .field("misses", &self.misses)
            .finish()
    }
}
