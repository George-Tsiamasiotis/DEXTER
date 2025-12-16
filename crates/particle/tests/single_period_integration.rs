use std::path::PathBuf;

use equilibrium::{
    bfields::*, currents::*, geometries::*, harmonics::*, perturbations::*, qfactors::*,
};
use particle::*;

#[test]
fn test_particle_single_period_integration() {
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
        theta0: 2.0,
        psip0: geometry.psip_wall() / 2.0,
        rho0: 1e-6,
        zeta0: 0.0,
        mu: 0.0,
    };

    // =================== Energy adaptive step

    let mut particle = Particle::new(&initial_conditions);
    let config = SinglePeriodConfig {
        method: SteppingMethod::EnergyAdaptiveStep,
        ..Default::default()
    };

    particle.single_period_integrate(&qfactor, &current, &bfield, &perturbation, &config);
    assert!(particle.status.is_single_period_integrated());
    assert!(particle.final_energy().is_finite());

    // =================== Error adaptive step

    let mut particle = Particle::new(&initial_conditions);
    let config = SinglePeriodConfig {
        method: SteppingMethod::ErrorAdaptiveStep,
        ..Default::default()
    };

    particle.single_period_integrate(&qfactor, &current, &bfield, &perturbation, &config);
    assert!(particle.status.is_single_period_integrated());
    assert!(particle.final_energy().is_finite());

    // =================== Few steps integration

    let mut particle = Particle::new(&initial_conditions);
    let config = SinglePeriodConfig {
        max_steps: 100,
        ..Default::default()
    };

    particle.single_period_integrate(&qfactor, &current, &bfield, &perturbation, &config);
    assert!(particle.status.is_timed_out());
    assert_eq!(particle.evolution.steps_taken(), 100);
    assert!(particle.final_energy().is_finite());

    let _ = format!("{:?}", &particle.evolution);
    particle.evolution.discard();
}
