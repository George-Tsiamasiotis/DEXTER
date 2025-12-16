use std::path::PathBuf;

use equilibrium::{
    bfields::*, currents::*, geometries::*, harmonics::*, perturbations::*, qfactors::*,
};
use ndarray::{Array1, Axis};
use particle::*;

#[test]
fn test_particle_integration() {
    let path = PathBuf::from("../equilibrium/lar_netcdf.nc");

    let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic")
        .build()
        .unwrap();
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = NcPerturbation::from_harmonics(&vec![
        NcHarmonicBuilder::new(&path, "steffen", 2, 1)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .build()
            .unwrap(),
    ]);

    let initial_conditions = InitialConditions {
        time0: 0.0,
        theta0: 0.0,
        psip0: geometry.psip_wall() / 2.0,
        rho0: 1e-4,
        zeta0: 0.0,
        mu: 0.0,
    };

    // =================== Energy adaptive step

    let mut particle = Particle::new(&initial_conditions);
    let config = IntegrationConfig {
        method: SteppingMethod::EnergyAdaptiveStep,
        ..Default::default()
    };

    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e6),
        &config,
    );
    assert!(particle.status.is_integrated());
    assert!(particle.final_energy().is_finite());

    // =================== Error adaptive step

    let mut particle = Particle::new(&initial_conditions);
    let config = IntegrationConfig {
        method: SteppingMethod::ErrorAdaptiveStep,
        ..Default::default()
    };

    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e6),
        &config,
    );
    assert!(particle.status.is_integrated());
    assert!(particle.final_energy().is_finite());

    // =================== Fixed step

    let mut particle = Particle::new(&initial_conditions);
    let stepsize = 1e-1;
    let config = IntegrationConfig {
        method: SteppingMethod::FixedStep(stepsize),
        ..Default::default()
    };

    let t_eval = (0.0, 10.0);
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        t_eval.clone(),
        &config,
    );
    assert_eq!(
        particle.evolution.steps_stored(),
        (t_eval.1 / stepsize) as usize + 1
    );
    let expected_diffs = Array1::from_elem(particle.evolution.steps_stored(), stepsize);
    let epsilon = stepsize * 1e-10;
    particle
        .evolution
        .time()
        .diff(1, Axis(0))
        .abs_diff_eq(&expected_diffs, epsilon);
    assert!(particle.status.is_integrated());
    assert!(particle.final_energy().is_finite());

    // =================== Few steps integration

    let mut particle = Particle::new(&initial_conditions);
    let config = IntegrationConfig {
        max_steps: 100,
        ..Default::default()
    };

    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e10),
        &config,
    );
    assert!(particle.status.is_timed_out());
    assert_eq!(particle.evolution.steps_taken(), 100);
    assert_eq!(particle.evolution.steps_stored(), 100);
    assert!(particle.final_energy().is_finite());

    let _ = format!("{:?}", &particle.evolution);
    particle.evolution.discard();
    let _ = format!("{:?}", particle);
}
