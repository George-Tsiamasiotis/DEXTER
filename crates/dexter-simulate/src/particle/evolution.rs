//! Stores the time evolution of a Particle.

use dexter_common::vec_to_array1D_getter_impl;
use ndarray::Array1;
use std::time::Duration;

use crate::state::GCState;

/// Time series for a Particle's orbit.
///
/// All values are in *Normalized units*.
#[non_exhaustive]
#[derive(Default)]
pub(crate) struct Evolution {
    pub(crate) t: Vec<f64>,
    pub(crate) psi: Vec<f64>,
    pub(crate) psip: Vec<f64>,
    pub(crate) theta: Vec<f64>,
    pub(crate) zeta: Vec<f64>,
    pub(crate) rho: Vec<f64>,
    pub(crate) mu: Vec<f64>,

    pub(crate) ptheta: Vec<f64>,
    pub(crate) pzeta: Vec<f64>,
    pub(crate) energy: Vec<f64>,

    pub(crate) duration: Duration,
    pub(crate) steps_taken: usize,

    /// Variance of the energy array.
    energy_var: Option<f64>,
}

impl Evolution {
    /// Resets all fields.
    pub(crate) fn reset(&self) -> Self {
        Self::default()
    }

    /// Pushes the calculated quantities of an *evaluated* [`GCState`] in the time series.
    pub(crate) fn push_state(&mut self, state: &GCState) {
        self.t.push(state.t);
        self.psi.push(state.psi);
        self.psip.push(state.psip);
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
        self.psi.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.mu.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.energy.shrink_to_fit();

        self.energy_var = Some(self.energy_array().var(1.0))
    }

    /// Discards the vecs, keeping all the other fields.
    #[allow(dead_code)]
    pub(crate) fn discard_arrays(&mut self) {
        *self = Self {
            duration: self.duration,
            steps_taken: self.steps_taken,
            energy_var: self.energy_var,
            ..Default::default()
        }
    }
}

/// Getters
impl Evolution {
    vec_to_array1D_getter_impl!(t_array, t, time);
    vec_to_array1D_getter_impl!(psi_array, psi, ψ);
    vec_to_array1D_getter_impl!(psip_array, psip, ψp);
    vec_to_array1D_getter_impl!(theta_array, theta, θ);
    vec_to_array1D_getter_impl!(zeta_array, zeta, ζ);
    vec_to_array1D_getter_impl!(rho_array, rho, ρ);
    vec_to_array1D_getter_impl!(mu_array, mu, μ);
    vec_to_array1D_getter_impl!(ptheta_array, ptheta, Pθ);
    vec_to_array1D_getter_impl!(pzeta_array, pzeta, Pζ);
    vec_to_array1D_getter_impl!(energy_array, energy, E);

    /// Returns the final time (last entry on the time array), if it exists.
    pub(crate) fn tf(&self) -> Option<f64> {
        self.t.last().copied()
    }

    pub(crate) fn steps_stored(&self) -> usize {
        self.t.len()
    }

    pub(crate) fn energy_var(&self) -> Option<f64> {
        self.energy_var
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
