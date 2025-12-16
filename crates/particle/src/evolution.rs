use std::fmt::Debug;
use std::time::Duration;

use common::array1D_getter_impl;
use ndarray::Array1;

use crate::State;

/// The initial capacity of the time series Vecs.
const EVOLUTION_INIT_CAPACITY: usize = 2000;

/// Time series for a Particle's orbit.
#[derive(Clone)]
pub struct Evolution {
    pub time: Vec<f64>,
    /// The `θ` angle time series.
    pub theta: Vec<f64>,
    /// The poloidal flux `ψp` time series.
    pub psip: Vec<f64>,
    /// The parallel gyroradius `ρ_{||}` time series.
    pub rho: Vec<f64>,
    /// The `ζ` angle time series.
    pub zeta: Vec<f64>,
    /// The toroidal flux `ψ` time series.
    pub psi: Vec<f64>,
    /// The canonical momentum `Pθ` time series.
    pub ptheta: Vec<f64>,
    /// The canonical momentum `Pζ` time series.
    pub pzeta: Vec<f64>,
    /// The energy time series.
    pub energy: Vec<f64>,
    /// The duration of the integration.
    pub duration: Duration,
    /// The total steps of the integration.
    pub steps_taken: usize,
    /// The steps stored in the time series.
    steps_stored: usize,
    /// Relative standard deviation of the energy time series (σ/μ).
    pub energy_std: f64,
}

impl Evolution {
    /// Creates an [`Evolution`], initializing the vecs with `capacity`.
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            time: Vec::with_capacity(capacity),
            theta: Vec::with_capacity(capacity),
            psip: Vec::with_capacity(capacity),
            rho: Vec::with_capacity(capacity),
            zeta: Vec::with_capacity(capacity),
            psi: Vec::with_capacity(capacity),
            ptheta: Vec::with_capacity(capacity),
            pzeta: Vec::with_capacity(capacity),
            energy: Vec::with_capacity(capacity),
            energy_std: f64::NAN,
            duration: Duration::default(),
            steps_taken: 0,
            steps_stored: 0,
        }
    }

    /// Returns the total steps of the integration.
    pub fn steps_taken(&self) -> usize {
        self.steps_taken
    }

    /// Returns the number of steps stored in each time series.
    pub fn steps_stored(&self) -> usize {
        self.steps_stored
    }

    /// Returns the final stored time.
    pub fn final_time(&self) -> Option<f64> {
        self.time.last().copied()
    }

    /// Pushes the variables of a [`State`] to the time series vecs.
    pub(crate) fn push_state(&mut self, state: &State) {
        self.time.push(state.time);
        self.theta.push(state.theta);
        self.psip.push(state.psip);
        self.rho.push(state.rho);
        self.zeta.push(state.zeta);
        self.psi.push(state.psi);
        self.ptheta.push(state.ptheta);
        self.pzeta.push(state.pzeta);
        self.energy.push(state.energy());
        self.steps_stored += 1;
    }

    /// Shrinks the vecs and calculates `energy_std`.
    pub(crate) fn finish(&mut self) {
        self.time.shrink_to_fit();
        self.theta.shrink_to_fit();
        self.psip.shrink_to_fit();
        self.rho.shrink_to_fit();
        self.zeta.shrink_to_fit();
        self.psi.shrink_to_fit();
        self.ptheta.shrink_to_fit();
        self.pzeta.shrink_to_fit();
        self.energy.shrink_to_fit();

        let energy_array = Array1::from_vec(self.energy.clone());
        self.energy_std = energy_array.std(0.0) / energy_array.mean().unwrap_or(f64::NAN);
    }

    /// Resets all arrays to the empty defaults, keeping all the other fields.
    ///
    /// Use this to free memory when dealing with many particles.
    pub fn discard(&mut self) {
        self.time = Vec::default();
        self.theta = Vec::default();
        self.psip = Vec::default();
        self.rho = Vec::default();
        self.zeta = Vec::default();
        self.psi = Vec::default();
        self.ptheta = Vec::default();
        self.pzeta = Vec::default();
        self.energy = Vec::default();
    }

    array1D_getter_impl!(time, time);
    array1D_getter_impl!(theta, theta);
    array1D_getter_impl!(psip, psip);
    array1D_getter_impl!(rho, rho);
    array1D_getter_impl!(zeta, zeta);
    array1D_getter_impl!(psi, psi);
    array1D_getter_impl!(ptheta, ptheta);
    array1D_getter_impl!(pzeta, pzeta);
    array1D_getter_impl!(energy, energy);
}

impl Debug for Evolution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Evolution")
            .field(
                "time",
                &format!(
                    "[{:.5}, {:.5}]",
                    self.time.first().unwrap_or(&f64::NAN),
                    self.time.last().unwrap_or(&f64::NAN),
                ),
            )
            .field("duration", &self.duration)
            .field("energy_std", &format!("{:.5}", self.energy_std))
            .field("steps taken", &self.steps_taken())
            .field("steps stored", &self.steps_stored())
            .finish()
    }
}

impl Default for Evolution {
    fn default() -> Self {
        Self::with_capacity(EVOLUTION_INIT_CAPACITY)
    }
}
