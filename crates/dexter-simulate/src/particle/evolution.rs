//! Stores the time evolution of a Particle.

use dexter_common::f64_slice_getter;
use std::time::Duration;

use crate::state::GCState;

/// Time series for a Particle's orbit.
///
/// All values are in *Normalized units*.
#[non_exhaustive]
#[derive(Default, Clone)]
pub(crate) struct Evolution {
    /// The time array.
    pub(crate) t: Vec<f64>,
    /// The `ψt` time-series.
    pub(crate) psit: Vec<f64>,
    /// The `ψp` time-series.
    pub(crate) psip: Vec<f64>,
    /// The `θ` time-series.
    pub(crate) theta: Vec<f64>,
    /// The `ζ` time-series.
    pub(crate) zeta: Vec<f64>,
    /// The `ρ` time-series.
    pub(crate) rho: Vec<f64>,
    /// The `μ` time-series. At the moment, `dμ/dt` is hardcoded to 0, but we treat it as a
    /// dynamical constant in case we add ξ-dependent perturbations.
    pub(crate) mu: Vec<f64>,

    /// The `Pθ` time-series.
    pub(crate) ptheta: Vec<f64>,
    /// The `Pζ` time-series.
    pub(crate) pzeta: Vec<f64>,
    /// The `E` time-series.
    pub(crate) energy: Vec<f64>,

    /// The total duration of the routine.
    pub(crate) duration: Duration,
    /// The total steps taken. Depending on the routine, this number might not be equal to the
    /// number of elements in the time arrays.
    pub(crate) steps_taken: usize,
}

impl Evolution {
    /// Resets all fields.
    pub(crate) fn reset(&mut self) {
        *self = Self::default()
    }

    /// Pushes the calculated quantities of an *evaluated* [`GCState`] in the time series.
    pub(crate) fn push_state(&mut self, state: &GCState) {
        self.t.push(state.t);
        self.psit.push(state.psit.value());
        self.psip.push(state.psip.value());
        self.theta.push(state.theta);
        self.zeta.push(state.zeta);
        self.rho.push(state.rho);
        self.mu.push(state.mu);
        self.ptheta.push(state.ptheta);
        self.pzeta.push(state.pzeta);
        self.energy.push(state.energy);
    }

    /// Shrinks the vecs and calculates the final energy variance.
    pub(crate) fn finish(&mut self) {
        self.t.shrink_to_fit();
        self.psit.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.mu.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.energy.shrink_to_fit();
    }

    /// Discards the vecs, keeping all the other fields.
    pub(crate) fn discard_arrays(&mut self) {
        let vectors = [
            &mut self.t,
            &mut self.psit,
            &mut self.psip,
            &mut self.theta,
            &mut self.zeta,
            &mut self.rho,
            &mut self.mu,
            &mut self.ptheta,
            &mut self.pzeta,
            &mut self.energy,
        ];
        for vector in vectors {
            *vector = Vec::new();
        }
    }
}

/// Getters.
impl Evolution {
    f64_slice_getter!(t_values, t, t);
    f64_slice_getter!(psit_values, psit, psit);
    f64_slice_getter!(psip_values, psip, psip);
    f64_slice_getter!(theta_values, theta, theta);
    f64_slice_getter!(zeta_values, zeta, zeta);
    f64_slice_getter!(rho_values, rho, rho);
    f64_slice_getter!(mu_values, mu, mu);
    f64_slice_getter!(ptheta_values, ptheta, ptheta);
    f64_slice_getter!(pzeta_values, pzeta, pzeta);
    f64_slice_getter!(energy_values, energy, energy);

    /// Returns the final time (last entry on the time array), if it exists.
    pub(crate) fn tf(&self) -> Option<f64> {
        self.t.last().copied()
    }

    /// Returns the number of steps stored in the time series.
    pub(crate) fn steps_stored(&self) -> usize {
        self.t.len()
    }
}

impl std::fmt::Debug for Evolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParticleEvolution")
            .field(
                "time",
                &format!(
                    "[{:.6}, {:.6}]",
                    self.t.first().unwrap_or(&f64::NAN),
                    self.t.last().unwrap_or(&f64::NAN),
                ),
            )
            .field("array length", &self.t.len())
            .field("duration", &self.duration)
            .field("steps_taken", &self.steps_taken)
            .field("steps_stored", &self.steps_stored())
            .finish()
    }
}
