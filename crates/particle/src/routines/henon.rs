//! Functions for calculating Poincare intersection States and period closures.

use std::f64::consts::{PI, TAU};

use equilibrium::{Bfield, Current, Perturbation, Qfactor};

use crate::Result;
use crate::{MappingParameters, PoincareSection};
use crate::{State, Stepper};

/// Creates a [`State`] on the modified system (6) from `state1`, which is a state of the normal
/// system. The modified state corresponds to the first state on the modified system.
pub(crate) fn calculate_mod_state1(state1: &State, section: &PoincareSection) -> State {
    // Do not evaluate the state!
    match section {
        PoincareSection::ConstTheta => {
            let kappa = 1.0 / state1.theta_dot;
            let dt_dtheta = kappa;
            let dpsip_dtheta = kappa * state1.psip_dot;
            let drho_dtheta = kappa * state1.rho_dot;
            let dzeta_dtheta = kappa * state1.zeta_dot;
            State {
                time: state1.theta,
                theta: state1.time,
                theta_dot: dt_dtheta,
                psip_dot: dpsip_dtheta,
                rho_dot: drho_dtheta,
                zeta_dot: dzeta_dtheta,
                hcache: state1.hcache.clone(),
                ..*state1
            }
        }
        PoincareSection::ConstZeta => {
            let kappa = 1.0 / state1.zeta_dot;
            let dtheta_dzeta = kappa * state1.theta_dot;
            let dpsip_dzeta = kappa * state1.psip_dot;
            let drho_dzeta = kappa * state1.rho_dot;
            let dt_dzeta = kappa;
            State {
                time: state1.zeta,
                zeta: state1.time,
                theta_dot: dtheta_dzeta,
                psip_dot: dpsip_dzeta,
                rho_dot: drho_dzeta,
                zeta_dot: dt_dzeta,
                hcache: state1.hcache.clone(),
                ..*state1
            }
        }
    }
}

/// Calculates the step size dτ that brings mod_state1 on the intersection surface.
pub(crate) fn calculate_mod_step(
    state1: &State,
    state2: &State,
    params: &MappingParameters,
) -> f64 {
    // WARN: This is needed to move the %2π pole when α happens to be a multiple of 2π. It seems to
    // work for most cases, but lets keep an eye on it.
    let pole = if params.alpha.abs() < 1e-2 { PI } else { 0.0 };
    match params.section {
        PoincareSection::ConstTheta => {
            let direction = (state2.theta - state1.theta).signum();
            direction * (params.alpha - (state1.theta + pole).rem_euclid(TAU) + pole)
        }
        PoincareSection::ConstZeta => {
            let direction = (state2.zeta - state1.zeta).signum();
            direction * (params.alpha - (state1.zeta + pole).rem_euclid(TAU) + pole)
        }
    }
}

/// Performs 1 step on the modified system (6) to calculate mod_state2, which sits exactly on the
/// intersection suface, **but corresponds to the modified system**.
pub(crate) fn calculate_mod_state2(
    qfactor: &impl Qfactor,
    current: &impl Current,
    bfield: &impl Bfield,
    perturbation: &impl Perturbation,
    mod_state1: State,
    dtau: f64,
) -> Result<State> {
    let mut mod_stepper = Stepper::new(&mod_state1);
    mod_stepper.start(dtau, qfactor, current, bfield, perturbation)?;
    {
        // NOTE: I am not sure why we must do this, but we must.
        // Probably equivalent to re-adjusting the step-size `h`, but only for the modified
        // system's independent variables, while keeping `h=dτ` for the 'τ' (modified system's
        // time) step size.
        mod_stepper.weights[0] = mod_state1.theta_dot;
        mod_stepper.weights[1] = mod_state1.psip_dot;
        mod_stepper.weights[2] = mod_state1.rho_dot;
        mod_stepper.weights[3] = mod_state1.zeta_dot;
    }
    let mod_state2 = mod_stepper.next_state(dtau);
    Ok(mod_state2)
}

/// Calculates the state of the original system exactly on the intersection surface, by converting
/// mod_state2 back on the original system.
pub(crate) fn calculate_intersection_state(
    qfactor: &impl Qfactor,
    current: &impl Current,
    bfield: &impl Bfield,
    perturbation: &impl Perturbation,
    params: &MappingParameters,
    mod_state2: State,
) -> Result<State> {
    match params.section {
        PoincareSection::ConstTheta => {
            let kappa = 1.0;
            let dt_dt = kappa;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dzeta_dt = kappa * mod_state2.zeta_dot;
            State {
                time: mod_state2.theta,
                theta: mod_state2.time,
                theta_dot: dt_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dzeta_dt,
                mu: mod_state2.mu,
                ..mod_state2
            }
        }
        PoincareSection::ConstZeta => {
            let kappa = 1.0;
            let dtheta_dt = kappa * mod_state2.theta_dot;
            let dpsip_dt = kappa * mod_state2.psip_dot;
            let drho_dt = kappa * mod_state2.rho_dot;
            let dt_dt = kappa;
            State {
                time: mod_state2.zeta,
                zeta: mod_state2.time,
                theta_dot: dtheta_dt,
                psip_dot: dpsip_dt,
                rho_dot: drho_dt,
                zeta_dot: dt_dt,
                mu: mod_state2.mu,
                ..mod_state2
            }
        }
    }
    .into_evaluated(qfactor, current, bfield, perturbation)
}

/// Checks when an angle has intersected with the surface at `angle`.
/// (source: seems to work)
pub(crate) fn intersected(old_angle: f64, new_angle: f64, surface_angle: f64) -> bool {
    let diff1 = new_angle - surface_angle;
    let diff2 = old_angle - surface_angle;
    // NOTE: Use `<=` here for the case `surface_angle == 0`, since the sine of angles very close
    // to 0 (but not very close to 2π, 4π, ...)  return exactly 0.0.
    ((diff1) / 2.0).sin() * ((diff2) / 2.0).sin() <= 0.0
}

#[cfg(test)]
mod test {
    use super::*;
    use std::f64::consts::{PI, TAU};

    #[test]
    fn test_intersected() {
        // No need to check for negative numbers, since the signs would negate each other.

        let eps = 1e-12;

        assert!(intersected(1.0 - eps, 1.0 + eps, 1.0));
        assert!(intersected(1.0f64.next_down(), 1.0f64.next_up(), 1.0));
        assert!(!intersected(1.0 - eps, 1.0 - 2.0 * eps, 1.0));
        assert!(!intersected(1.0 + eps, 1.0 + 2.0 * eps, 1.0));
        assert!(!intersected(1.0 - eps, 1.0 - eps, 1.0));
        assert!(!intersected(1.0 + eps, 1.0 + eps, 1.0));
        assert!(!intersected(1.0f64.next_down(), 1.0f64.next_down(), 1.0));
        assert!(!intersected(1.0f64.next_up(), 1.0f64.next_up(), 1.0));

        assert!(intersected(10.0 - eps, 10.0 + eps, 10.0));
        assert!(intersected(10.0f64.next_down(), 10.0f64.next_up(), 10.0));
        assert!(!intersected(10.0 - eps, 10.0 - 2.0 * eps, 10.0));
        assert!(!intersected(10.0 + eps, 10.0 + 2.0 * eps, 10.0));
        assert!(!intersected(10.0 - eps, 10.0 - eps, 10.0));
        assert!(!intersected(10.0 + eps, 10.0 + eps, 10.0));
        assert!(!intersected(10.0f64.next_down(), 10.0f64.next_down(), 10.0));
        assert!(!intersected(10.0f64.next_up(), 10.0f64.next_up(), 10.0));

        assert!(intersected(PI - eps, PI + eps, PI));
        assert!(intersected(PI.next_down(), PI.next_up(), PI));
        assert!(!intersected(PI - eps, PI - 2.0 * eps, PI));
        assert!(!intersected(PI + eps, PI + 2.0 * eps, PI));
        assert!(!intersected(PI - eps, PI - eps, PI));
        assert!(!intersected(PI + eps, PI + eps, PI));
        assert!(!intersected(PI.next_down(), PI.next_down(), PI));
        assert!(!intersected(PI.next_up(), PI.next_up(), PI));

        assert!(intersected(TAU - eps, TAU + eps, TAU));
        assert!(intersected(TAU.next_down(), TAU.next_up(), TAU));
        assert!(!intersected(TAU - eps, TAU - 2.0 * eps, TAU));
        assert!(!intersected(TAU + eps, TAU + 2.0 * eps, TAU));
        assert!(!intersected(TAU - eps, TAU - eps, TAU));
        assert!(!intersected(TAU + eps, TAU + eps, TAU));
        assert!(!intersected(TAU.next_down(), TAU.next_down(), TAU));
        assert!(!intersected(TAU.next_up(), TAU.next_up(), TAU));

        assert!(intersected(0.0 - eps, 0.0 + eps, 0.0));
        assert!(intersected(0.0f64.next_down(), 0.0f64.next_up(), 0.0));
        assert!(!intersected(0.0 - eps, 0.0 - 2.0 * eps, 0.0));
        assert!(!intersected(0.0 + eps, 0.0 + 2.0 * eps, 0.0));
        assert!(!intersected(0.0 - eps, 0.0 - eps, 0.0));
        assert!(!intersected(0.0 + eps, 0.0 + eps, 0.0));
        // VERY corner case since all arguements have a sine of 0.0. Not worth working around. If
        // a particle is unlucky enough to reach this point, it would just be rejected from the
        // accuracy test.
        // assert!(!intersected(0.0f64.next_down(), 0.0f64.next_down(), 0.0));
        // assert!(!intersected(0.0f64.next_up(), 0.0f64.next_up(), 0.0));

        assert!(intersected(
            (2.0 * PI + PI).next_down(),
            (2.0 * PI + PI).next_up(),
            PI
        ));

        assert!(intersected(TAU.next_down(), TAU.next_up(), TAU));
        assert!(intersected(
            (2.0 * PI + TAU).next_down(),
            (2.0 * PI + TAU).next_up(),
            TAU
        ));

        assert!(!intersected(TAU - eps, TAU + eps, PI));
        assert!(!intersected(PI - eps, PI + eps, TAU));
        assert!(!intersected(PI - eps, PI + eps, PI / 2.0));
        assert!(!intersected(PI / 2.0 - eps, PI / 2.0 + eps, TAU));
    }
}
