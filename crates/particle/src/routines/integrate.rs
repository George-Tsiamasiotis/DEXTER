//! Integration of a [`Particle`] for a specific time interval.

use std::time::Instant;

use crate::Time;
use crate::{Evolution, Particle, State, Stepper};
use crate::{ParticleError, Result};

use config::*;
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};

/// Integrates the particle for a specific time interval.
pub(crate) fn integrate(
    particle: &mut Particle,
    qfactor: &Qfactor,
    bfield: &Bfield,
    currents: &Currents,
    perturbation: &Perturbation,
    t_eval: (Time, Time),
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

    // ==================== Main loop

    loop {
        if particle.evolution.steps_taken() == MAX_STEPS {
            return Err(ParticleError::TimedOut(start.elapsed()));
        }
        if particle.evolution.final_time().unwrap_or(t_eval.0) > t_eval.1 {
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        stepper.start(dt, qfactor, bfield, currents, perturbation)?;
        dt = stepper.calculate_optimal_step(dt)?;
        state2 = stepper
            .next_state(dt)
            .into_evaluated(qfactor, currents, bfield, perturbation)?;
        particle.evolution.steps_taken += 1;

        // Store and continue
        particle.evolution.push_state(&state2);
        state1 = state2;
    }

    // ==================== Finalization

    particle.final_state = state1.into_evaluated(qfactor, currents, bfield, perturbation)?;
    particle.evolution.finish();
    particle.calculate_orbit_type();
    particle.evolution.duration = start.elapsed();
    Ok(())
}
