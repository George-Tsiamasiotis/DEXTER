//! Integration of a [`Particle`] for a single period and calculation of ωθ, ωζ and qkinetic.

use std::f64::consts::PI;
use std::time::Instant;

use crate::routines::henon::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};
use crate::{Evolution, MappingParameters, Particle, PoincareSection, State, Stepper};
use crate::{ParticleError, Result};
use config::{RKF45_FIRST_STEP, SINGLE_PERIOD_MAX_STEPS};

use equilibrium::{Bfield, Currents, Perturbation, Qfactor};

/// A particle's `ωθ`, `ωζ` and `qkinetic`.
#[derive(Default, Clone)]
pub struct Frequencies {
    omega_theta: Option<f64>,
    omega_zeta: Option<f64>,
    qkinetic: Option<f64>,
    pub(crate) psip_intersections: (usize, Vec<usize>),
    pub(crate) theta_intersections: (usize, Vec<usize>),
}

impl Frequencies {
    /// Updates ωθ to a new value, and updates qkinetic accordingly.
    pub(crate) fn update_omega_theta(&mut self, value: f64) {
        self.omega_theta = Some(value);
        self.update_qkinetic();
    }

    /// Updates ωθ to a new value, and updates qkinetic accordingly.
    pub(crate) fn update_omega_zeta(&mut self, value: f64) {
        self.omega_zeta = Some(value);
        self.update_qkinetic();
    }

    /// Sets qkinetic to ωζ/ωθ if both fields are Some(), otherwise None.
    fn update_qkinetic(&mut self) {
        self.qkinetic = if self.omega_theta.is_some() && self.omega_zeta.is_some() {
            Some(
                self.omega_zeta.expect("already checked")
                    / self.omega_theta.expect("already checked"),
            )
        } else {
            None
        }
    }

    pub(crate) fn update_theta_intersections(&mut self, step_num: usize) {
        self.theta_intersections.0 += 1;
        self.theta_intersections.1.push(step_num);
    }

    pub(crate) fn update_psip_intersections(&mut self, step_num: usize) {
        self.psip_intersections.0 += 1;
        self.psip_intersections.1.push(step_num);
    }

    pub fn omega_theta(&self) -> Option<f64> {
        self.omega_theta
    }

    pub fn omega_zeta(&self) -> Option<f64> {
        self.omega_zeta
    }

    pub fn qkinetic(&self) -> Option<f64> {
        self.qkinetic
    }
}

// ===============================================================================================

/// Integrates the particle for 1 `θ-ψp` period.
pub(crate) fn close_theta_period(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    currents: &Currents,
    perturbation: &Perturbation,
) -> Result<()> {
    // ==================== Setup

    let start = Instant::now();
    particle.evolution = Evolution::default();
    particle
        .initial_state
        .evaluate(qfactor, currents, bfield, perturbation)?;
    let mut state1 = particle.initial_state.clone();
    let mut state2: State;
    let mut dt = RKF45_FIRST_STEP;

    let theta0 = particle.initial_state.theta;
    let psip0 = particle.initial_state.psip;
    let zeta0 = particle.initial_state.zeta;
    let time0 = particle.initial_state.time;
    let theta_dot0 = particle.initial_state.theta_dot;
    let psip_dot0 = particle.initial_state.psip_dot;

    // ==================== Main loop

    loop {
        if particle.evolution.steps_taken() == SINGLE_PERIOD_MAX_STEPS {
            return Err(ParticleError::TimedOut(start.elapsed()));
        }

        // Perform a step.
        let mut stepper = Stepper::new(&state1);
        stepper.start(dt, qfactor, bfield, currents, perturbation)?;
        dt = stepper.calculate_optimal_step(dt)?;
        state2 = stepper
            .next_state(dt)
            .into_evaluated(qfactor, currents, bfield, perturbation)?;
        particle.evolution.steps_taken += 1;

        // Prevent particle from stopping at the first step
        if particle.evolution.steps_stored() < 2 {
            state1 = state2;
            continue;
        }

        // Check if close to a period.
        let old_theta = state1.theta;
        let new_theta = state2.theta;
        let old_psip = state1.psip;
        let new_psip = state2.psip;
        let old_theta_dot = state1.theta_dot;
        let old_psip_dot = state1.psip_dot;

        let mut psip_intersected = false;
        let mut theta_intersected = false;
        let mut same_directions = false;

        // Check for intersections
        // TODO: An additional proximity check is needed for ψp, for the case where 𝜃₀=0,π, where
        // ψp tends to be a local extremum.
        if intersected(old_psip, new_psip, psip0) {
            particle
                .frequencies
                .update_psip_intersections(particle.evolution.steps_taken());
            psip_intersected = true;
        };
        if intersected(old_theta, new_theta, theta0) {
            particle
                .frequencies
                .update_theta_intersections(particle.evolution.steps_taken());
            theta_intersected = true;
        };

        // Check for direction
        if (old_theta_dot.signum() == theta_dot0.signum())
            && (old_psip_dot.signum() == psip_dot0.signum())
        {
            same_directions = true
        }

        if theta_intersected && psip_intersected && same_directions {
            // Hénon's trick.
            // If the particle intersected the `θ0-ψp0` point, go back to `state1` and find
            // the exact time step needed to complete the period. This step brings `θ` to
            // its initial value *exactly*, but not `ψp`, although the difference is negligible
            let params = MappingParameters {
                section: PoincareSection::ConstTheta,
                alpha: theta0,
                intersections: 1,
            };
            let mod_state1 = calculate_mod_state1(&state1, &params.section);
            let dtau = calculate_mod_step(&state1, &state2, &params);
            let mod_state2 =
                calculate_mod_state2(qfactor, bfield, currents, perturbation, mod_state1, dtau)?;
            let intersection_state = calculate_intersection_state(
                qfactor,
                bfield,
                currents,
                perturbation,
                &params,
                mod_state2,
            )?;

            // ================ Frequencies calculation

            // NOTE:
            // >>> The orbit frequency ωζ corresponds to the bounce/transit averaged rate of
            // >>> toroidal precession Δζ/Tω.
            let omega_period = intersection_state.time - time0;
            let omega_theta = 2.0 * PI / omega_period;
            let dzeta = intersection_state.zeta - zeta0;
            let omega_zeta = dzeta / omega_period;

            particle.frequencies.update_omega_theta(omega_theta);
            particle.frequencies.update_omega_zeta(omega_zeta);
            particle.evolution.push_state(&intersection_state);
            particle.final_state = intersection_state;
            break;
        }
        // Keep going if not close to a period.
        state1 = state2;
    }

    // ==================== Finalization

    particle.final_state = state1.into_evaluated(qfactor, currents, bfield, perturbation)?;
    particle.evolution.finish();
    particle.calculate_orbit_type();
    particle.evolution.duration = start.elapsed();
    Ok(())
}

// ===============================================================================================

impl std::fmt::Debug for Frequencies {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn stringify(o: Option<f64>) -> String {
            if let Some(value) = o {
                format!("{:.7}", value)
            } else {
                String::from("Not calculated")
            }
        }
        let theta_intersections = format!("{:?}", &self.theta_intersections);
        let psip_intersections = format!("{:?}", &self.psip_intersections);

        f.debug_struct("Frequencies")
            .field("omega_theta", &stringify(self.omega_theta))
            .field("omega_zeta", &stringify(self.omega_zeta))
            .field("qkinetic", &stringify(self.qkinetic))
            .field("ψp-intersections", &psip_intersections)
            .field("θ-intersections", &theta_intersections)
            .finish()
    }
}
