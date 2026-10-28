//! Common testing utilities.

use dexter_simulate::*;

/// Checks that the time series calculated from an integration routine are valid.
#[allow(dead_code, reason = "used in tests")]
pub(crate) fn check_integrated_particle_time_series(particle: &Particle) {
    let _: &InitialConditions = particle.initial_conditions();
    assert_eq!(particle.t_values().len(), particle.steps_stored());
    assert_eq!(particle.psit_values().len(), particle.steps_stored());
    assert_eq!(particle.psip_values().len(), particle.steps_stored());
    assert_eq!(particle.theta_values().len(), particle.steps_stored());
    assert_eq!(particle.zeta_values().len(), particle.steps_stored());
    assert_eq!(particle.rho_values().len(), particle.steps_stored());
    assert_eq!(particle.mu_values().len(), particle.steps_stored());
    assert_eq!(particle.ptheta_values().len(), particle.steps_stored());
    assert_eq!(particle.pzeta_values().len(), particle.steps_stored());
    assert_eq!(particle.energy_values().len(), particle.steps_stored());

    assert!(particle.t_values().iter().all(|x| x.is_finite()));
    assert!(particle.psit_values().iter().all(|x| x.is_finite()));
    assert!(particle.psip_values().iter().all(|x| x.is_finite()));
    assert!(particle.theta_values().iter().all(|x| x.is_finite()));
    assert!(particle.zeta_values().iter().all(|x| x.is_finite()));
    assert!(particle.rho_values().iter().all(|x| x.is_finite()));
    assert!(particle.mu_values().iter().all(|x| x.is_finite()));
    assert!(particle.ptheta_values().iter().all(|x| x.is_finite()));
    assert!(particle.pzeta_values().iter().all(|x| x.is_finite()));
    assert!(particle.energy_values().iter().all(|x| x.is_finite()));
}
