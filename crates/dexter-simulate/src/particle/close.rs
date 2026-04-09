//! Integration of a [`Particle`] for a specific amound of full `θ-ψ` periods.

use std::f64::consts::TAU;
use std::time::Instant;

use approx::relative_eq;
use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, HarmonicCache, Qfactor};
use ndarray_stats::QuantileExt;

use crate::particle::intersect::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};
use crate::particle::{EqObjects, IntegrationCaches, Particle, ParticleCacheStats};
use crate::solve::{SolverParams, Stepper};
use crate::state::GCState;
use crate::{Frequencies, IntersectParams, OrbitType};

use super::{IntegrationStatus, Intersection};

/// The relative tolerance of the `ψ-ψ0` check.
const PSI_RELATIVE_TOLERANCE: f64 = 1e-6;
/// See [`OrbitType::TrappedStagnated`].
const TRAPPED_STAGNATED_CLASSIFICATION_CUTOFF: f64 = 0.999;

// ===============================================================================================

/// We dont want this function to return an error; Instead, we want to set a corresponding
/// [`IntegrationStatus`] variant for each possible error.
pub(super) fn close<Q, C, B, H>(
    particle: &mut Particle,
    objects: &EqObjects<Q, C, B, H>,
    periods: usize,
    solver_params: &SolverParams,
) where
    Q: Qfactor + FluxCommute,
    C: Current,
    B: Bfield,
    H: Harmonic,
{
    // =============== Setup

    let start = Instant::now();
    particle.evolution.reset();
    let mut caches = IntegrationCaches::<H::Cache> {
        harmonic_caches: objects.perturbation.generate_caches(),
        ..Default::default()
    };
    let mut mod_caches = IntegrationCaches::<H::Cache> {
        harmonic_caches: objects.perturbation.generate_caches(),
        ..Default::default()
    };

    // Return early if the initial flux happens to be exactly 0.0 or out of bounds.
    if particle.initial_conditions().flux0.value() == 0.0 {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    }
    if particle.initial_conditions.finalize(objects).is_err() {
        particle.integration_status = IntegrationStatus::InvalidInitialConditions;
        return;
    }
    let Ok(mut state1) = GCState::new(&particle.initial_conditions, objects, &mut caches) else {
        particle.integration_status = IntegrationStatus::OutOfBoundsInitialization;
        return;
    };

    // `state1` has been evaluated on the initial point
    particle.initial_energy = Some(state1.energy());
    let theta_dot0_sign = state1.theta_dot.signum();
    let theta0 = state1.theta;
    let psi0 = state1.psi;

    let mut state2: GCState;
    let mut dt = solver_params.first_step;
    let mut closed_periods: usize = 0;

    // Phony intersection parameters to be used in the modified system. The `turns` field is
    // ignored.
    let intersect_params = &IntersectParams::new(Intersection::ConstTheta, theta0, 0);

    // =============== Main loop

    loop {
        if particle.evolution.steps_taken == solver_params.max_steps {
            particle.integration_status = match closed_periods {
                0 => IntegrationStatus::TimedOut(start.elapsed()),
                1.. => IntegrationStatus::ClosedPeriods(closed_periods),
            };
            break;
        }

        // Perform a step
        let mut stepper = Stepper::new(&state1);
        state2 = if let Ok(state) = stepper
            .start(dt, objects, &mut caches)
            .inspect(|_| dt = stepper.calculate_optimal_step(dt, solver_params))
            .and_then(|_| stepper.next_state(dt, objects, &mut caches))
        {
            state
        } else {
            // `start()` and `next_state()` can only fail if an evaluation is out of bounds.
            particle.integration_status = IntegrationStatus::Escaped;
            break;
        };

        // Avoid stopping the particle at the start of the integration. When the `stepping_method`
        // is `EnergyAdaptiveStep`, the step size grows almost exponentially, so 10 steps should be
        // enough for `ψ` to move sufficiently far from `ψ0`.
        if closed_period(&state1, &state2, theta0, psi0, theta_dot0_sign)
            && particle.steps_taken() > 10
        {
            closed_periods += 1;
        };

        if closed_periods == periods {
            particle.integration_status = IntegrationStatus::ClosedPeriods(periods);
            particle.evolution.push_state(&state1);
            break;
        }

        // Store and continue
        particle.evolution.push_state(&state1);
        particle.evolution.steps_taken += 1;
        state1 = state2;
    }

    // =============== Hénon trick to close the period

    // When we are about to close the final period, we use Hénon's trick to take one more step
    // and land exactly on the initial point. This significantly improves the accuracy of the
    // calculations of the frequencies.
    if particle.integration_status == IntegrationStatus::ClosedPeriods(periods) {
        // Switch to the modified system
        let mod_state1 = calculate_mod_state1(&state1, &intersect_params.intersection);

        // Time step to accurately close the orbit
        let dtau = calculate_mod_step(&state1, intersect_params);

        // Perform the step on the modified system
        // Can only fail if an evaluation is out of bounds
        //
        // If `mod_state2` was calculated correctly, switch back to the normal system, calculate
        // `intersection_state` and store it
        match calculate_mod_state2(objects, &mod_state1, dtau, &mut mod_caches) {
            Ok(mod_state2) => {
                match calculate_intersection_state(
                    objects,
                    &mod_state2,
                    intersect_params,
                    &mut mod_caches,
                ) {
                    Ok(intersection_state) => {
                        particle.evolution.push_state(&intersection_state);
                        particle.integration_status =
                            IntegrationStatus::ClosedPeriods(closed_periods)
                    }
                    Err(_) => {
                        particle.integration_status = IntegrationStatus::ModStateEscaped;
                    }
                }
            }
            Err(_) => {
                particle.integration_status = IntegrationStatus::ModStateEscaped;
            }
        };
    }

    // =============== Finalize

    classify_orbit(particle);
    calculate_frequencies(particle, closed_periods);
    particle.evolution.duration = start.elapsed();
    particle.final_energy = Some(state1.energy());
    particle.evolution.finish();
    particle.stats = ParticleCacheStats {
        psi_acc: caches.psi_acc,
        psip_acc: caches.psip_acc,
        theta_acc: caches.theta_acc,
        harmonic_cache_hits: caches.harmonic_caches.iter().map(H::Cache::hits).sum(),
        harmonic_cache_misses: caches.harmonic_caches.iter().map(H::Cache::misses).sum(),
    };
}

/// Checks if `state1` and `state2` stand 'left and right' of the period closing point. This
/// indicates that the particle reached its starting point.
#[expect(clippy::float_cmp, reason = "sign check")]
fn closed_period(
    state1: &GCState,
    state2: &GCState,
    theta0: f64,
    psi0: f64,
    theta_dot0_sign: f64,
) -> bool {
    // Short-circuit the `ψ` check since `intersected()` can be significantly more expensive to
    // calculate.
    //
    // It should probably take into account the `θ` direction too.
    (relative_eq!(state1.psi, psi0, epsilon = PSI_RELATIVE_TOLERANCE)
        || intersected(state1.psi, state2.psi, psi0))
        && intersected(state1.theta, state2.theta, theta0)
        && state1.theta_dot.signum() == theta_dot0_sign
}

/// Classifies the particle's [`OrbitType`].
///
/// See [`OrbitType`]'s variants' documentation for the definitions.
fn classify_orbit(particle: &mut Particle) {
    // `min()` and `max()` can only fail if `theta_array` contains `NaN`s
    let Ok(theta_min) = particle.theta_array().min().copied() else {
        particle.orbit_type = OrbitType::Unclassified;
        return;
    };
    let Ok(theta_max) = particle.theta_array().max().copied() else {
        particle.orbit_type = OrbitType::Unclassified;
        return;
    };
    if (theta_max - theta_min).abs() < TAU * TRAPPED_STAGNATED_CLASSIFICATION_CUTOFF {
        particle.orbit_type = OrbitType::TrappedStagnated;
        return;
    }

    if particle
        .rho_array()
        .iter()
        .all(|rho| rho.is_sign_positive())
    {
        particle.orbit_type = OrbitType::CoPassing
    } else if particle
        .rho_array()
        .iter()
        .all(|rho| rho.is_sign_negative())
    {
        particle.orbit_type = OrbitType::CuPassing
    } else {
        particle.orbit_type = OrbitType::Unclassified
    }
}

/// Calculates the final `ωθ`, `ωζ` and `qkinetic` and stores them in the particle's
/// `frequencies` field.
///
/// Also takes account the number of periods closed.
fn calculate_frequencies(particle: &mut Particle, periods: usize) {
    let periodsf64 = periods as f64;
    // Use `NaN`for particles with invalid initial conditions or out of bounds initialization.
    let t0 = particle.t_array().first().copied().unwrap_or(f64::NAN);
    let tf = particle.t_array().last().copied().unwrap_or(f64::NAN);
    let theta_period = tf - t0;

    let zeta0 = particle.zeta_array().first().copied().unwrap_or(f64::NAN);
    let zetaf = particle.zeta_array().last().copied().unwrap_or(f64::NAN);
    let dzeta = zetaf - zeta0;

    let omega_theta = TAU / theta_period / periodsf64;
    let omega_zeta = dzeta / theta_period / periodsf64;
    let qkinetic = omega_zeta / omega_theta;

    particle.frequencies = Frequencies {
        omega_theta: Some(omega_theta),
        omega_zeta: Some(omega_zeta),
        qkinetic: Some(qkinetic),
    }
}
