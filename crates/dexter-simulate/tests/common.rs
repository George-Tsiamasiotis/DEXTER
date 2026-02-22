use dexter_simulate::Particle;

pub fn check_integrated_particle_arrays(particle: &Particle) {
    particle.initial_conditions();
    assert_eq!(particle.t_array().len(), particle.steps_stored());
    assert_eq!(particle.psi_array().len(), particle.steps_stored());
    assert_eq!(particle.psip_array().len(), particle.steps_stored());
    assert_eq!(particle.theta_array().len(), particle.steps_stored());
    assert_eq!(particle.zeta_array().len(), particle.steps_stored());
    assert_eq!(particle.rho_array().len(), particle.steps_stored());
    assert_eq!(particle.mu_array().len(), particle.steps_stored());
    assert_eq!(particle.ptheta_array().len(), particle.steps_stored());
    assert_eq!(particle.pzeta_array().len(), particle.steps_stored());
    assert_eq!(particle.energy_array().len(), particle.steps_stored());

    assert!(particle.t_array().iter().all(|x| x.is_finite()));
    assert!(particle.psi_array().iter().all(|x| x.is_finite()));
    assert!(particle.psip_array().iter().all(|x| x.is_finite()));
    assert!(particle.theta_array().iter().all(|x| x.is_finite()));
    assert!(particle.zeta_array().iter().all(|x| x.is_finite()));
    assert!(particle.rho_array().iter().all(|x| x.is_finite()));
    assert!(particle.mu_array().iter().all(|x| x.is_finite()));
    assert!(particle.ptheta_array().iter().all(|x| x.is_finite()));
    assert!(particle.pzeta_array().iter().all(|x| x.is_finite()));
    assert!(particle.energy_array().iter().all(|x| x.is_finite()));
}
