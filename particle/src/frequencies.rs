use std::f64::consts::PI;

use config::{MAX_STEPS, RKF45_FIRST_STEP};
use equilibrium::{Bfield, Currents, Perturbation, Qfactor};
use safe_unwrap::safe_unwrap;

use crate::Omega::{OmegaTheta, OmegaZeta};
use crate::PoincareSection;
use crate::Result;
use crate::{MappingParameters, Particle, State, Stepper};

use crate::henon::{
    calculate_intersection_state, calculate_mod_state1, calculate_mod_state2, calculate_mod_step,
    intersected,
};

/// Frequency selector.
pub(crate) enum Omega {
    OmegaTheta,
    #[allow(unused)]
    OmegaZeta,
}

/// A particle's `ωθ`, `ωζ` and `qkinetic`.
#[derive(Default, Clone)]
pub struct Frequencies {
    omega_theta: Option<f64>,
    omega_zeta: Option<f64>,
    qkinetic: Option<f64>,
}

impl Frequencies {
    /// Updates a frequency to a new value, and updates qkinetic accordingly.
    pub(crate) fn update(&mut self, which_omega: &Omega, value: f64) {
        use Omega::*;
        match which_omega {
            OmegaTheta => self.omega_theta = Some(value),
            OmegaZeta => self.omega_zeta = Some(value),
        };
        self.update_qkinetic();
    }

    /// Sets qkinetic to ωζ/ωθ if both fields are Some(), otherwise None.
    fn update_qkinetic(&mut self) {
        self.qkinetic = if self.omega_theta.is_some() && self.omega_zeta.is_some() {
            Some(
                safe_unwrap!("already checked", self.omega_zeta)
                    / safe_unwrap!("already checked", self.omega_theta),
            )
        } else {
            None
        }
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
    // Last two states of the particle.
    let mut state1 = particle.initial_state.clone(); // Already evaluated
    let mut state2: State;

    let theta0 = particle.initial_state.theta;
    let psip0 = particle.initial_state.psip;
    let zeta0 = particle.initial_state.zeta;
    let time0 = particle.initial_state.time;

    let mut dt = RKF45_FIRST_STEP;
    while particle.evolution.steps_taken() <= MAX_STEPS {
        particle.evolution.push_state(&state1);
        particle.evolution.steps_taken += 1;

        // Perform a step
        let mut stepper = Stepper::default();
        stepper.init(&state1);
        stepper.start(dt, qfactor, bfield, currents, perturbation)?;
        dt = stepper.calculate_optimal_step(dt)?;
        state2 = stepper.next_state(dt);

        state2.evaluate(qfactor, currents, bfield, perturbation)?;

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
        // Use `intersected` rather than `is_close` checks to avoid stopping the particles
        // immediately and hardcoding tolerances
        // Intersected() for the flux might be unnecessary here, but its safe.
        if intersected(old_psip, new_psip, psip0) && intersected(old_theta, new_theta, theta0) {
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
            )?
            .into_evaluated(qfactor, currents, bfield, perturbation)?;

            #[allow(non_snake_case)]
            let Tomega = intersection_state.time - time0;
            let omega_theta = 2.0 * PI / Tomega;

            // NOTE:
            // >>> The orbit frequency ωζ corresponds to the bounce/transit averaged rate of
            // >>> toroidal precession Δζ/Tω.
            let dzeta = intersection_state.zeta - zeta0;
            let omega_zeta = dzeta / Tomega;

            particle.frequencies.update(&OmegaTheta, omega_theta);
            particle.frequencies.update(&OmegaZeta, omega_zeta);
            particle.evolution.push_state(&intersection_state);
            particle.final_state = intersection_state;
            particle.status = crate::IntegrationStatus::SinglePeriodIntegrated;
            break;
        }
        // Keep going if not close to a period.
        state1 = state2;
    }

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

        f.debug_struct("Frequencies")
            .field("omega_theta", &stringify(self.omega_theta))
            .field("omega_zeta", &stringify(self.omega_zeta))
            .field("qkinetic", &stringify(self.qkinetic))
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use Omega::*;

    #[test]
    fn test_frequencies_update_default() {
        let mut freq = Frequencies::default();

        assert!(freq.omega_theta.is_none());
        assert!(freq.omega_zeta.is_none());
        assert!(freq.qkinetic.is_none());

        freq.update(&OmegaTheta, 2.0);

        assert_eq!(freq.omega_theta.unwrap(), 2.0);
        assert!(freq.omega_zeta.is_none());
        assert!(freq.qkinetic.is_none());

        freq.update(&OmegaZeta, 20.0);

        assert_eq!(freq.omega_theta.unwrap(), 2.0);
        assert_eq!(freq.omega_zeta.unwrap(), 20.0);
        assert_eq!(freq.qkinetic.unwrap(), 10.0);

        freq.update(&OmegaTheta, 4.0);

        assert_eq!(freq.omega_theta.unwrap(), 4.0);
        assert_eq!(freq.omega_zeta.unwrap(), 20.0);
        assert_eq!(freq.qkinetic.unwrap(), 5.0);

        freq.update(&OmegaZeta, 100.0);

        assert_eq!(freq.omega_theta.unwrap(), 4.0);
        assert_eq!(freq.omega_zeta.unwrap(), 100.0);
        assert_eq!(freq.qkinetic.unwrap(), 25.0);
    }
}
