//! Tests corner cases for unwanted particle integration results

mod common;
use dexter_machine::*;
use dexter_simulate::*;

#[test]
fn time_out() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.4));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.1), 0.0, 0.0, 1e-4, 1e-6);
    let solver_params = SolverParams {
        max_steps: 10,
        ..Default::default()
    };

    // Integrate
    let mut particle = Particle::new(&initial);
    particle.integrate(machine, (0.0, 1e5), &solver_params);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));

    // Intersect
    let mut particle = Particle::new(&initial);
    let intersect_params =
        IntersectParams::new(Intersection::ConstTheta, 1.0, 10, Directionality::default());
    particle.intersect(machine, &intersect_params, &solver_params);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));

    // Close
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &solver_params);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::TimedOut(..)
    ));
}

#[test]
fn out_of_bounds_initialization() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.4));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(1.0), 0.0, 0.0, 1e-4, 1e-6);

    // Integrate
    let mut particle = Particle::new(&initial);
    particle.integrate(machine, (0.0, 1e5), &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("OutOfBoundsInitialization".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::OutOfBoundsInitialization
    );

    // Intersect
    let mut particle = Particle::new(&initial);
    let intersect_params =
        IntersectParams::new(Intersection::ConstTheta, 0.0, 10, Directionality::default());
    particle.intersect(machine, &intersect_params, &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("OutOfBoundsInitialization".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::OutOfBoundsInitialization
    );

    // Close
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("OutOfBoundsInitialization".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::OutOfBoundsInitialization
    );

    // Classify
    let mut particle = Particle::new(&initial);
    particle.classify(machine);
    assert_eq!(particle.integration_status(), IntegrationStatus::None);
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::OutOfBoundsInitialization
    );
}

#[test]
fn bad_magnetic_flux_initialization() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.4));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial =
        InitialConditions::boozer(0.0, MagneticFlux::Poloidal(0.02), 0.0, 0.0, 1e-4, 1e-6);

    // Integrate
    let mut particle = Particle::new(&initial);
    particle.integrate(machine, (0.0, 1e5), &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("BadMagneticFlux".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::BadMagneticFlux
    );

    // Intersect
    let mut particle = Particle::new(&initial);
    let intersect_params =
        IntersectParams::new(Intersection::ConstTheta, 0.0, 10, Directionality::default());
    particle.intersect(machine, &intersect_params, &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("BadMagneticFlux".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::BadMagneticFlux
    );

    // Close
    let mut particle = Particle::new(&initial);
    particle.close(machine, 1, &SolverParams::default());
    assert_eq!(
        particle.integration_status(),
        IntegrationStatus::Failed("BadMagneticFlux".into())
    );
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::BadMagneticFlux
    );

    // Classify
    let mut particle = Particle::new(&initial);
    particle.classify(machine);
    assert_eq!(particle.integration_status(), IntegrationStatus::None);
    assert_eq!(
        particle.initial_conditions().status(),
        InitialConditionsStatus::BadMagneticFlux
    );
}

#[test]
fn intersected_time_out() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.5));
    let current = LarCurrent;
    let bfield = LarBfield;
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let initial = InitialConditions::boozer(0.0, MagneticFlux::Toroidal(0.1), 0.0, 0.0, 1e-4, 0.0);
    let solver_params = SolverParams {
        max_steps: 10000,
        ..Default::default()
    };

    let mut particle = Particle::new(&initial);
    let intersect_params = IntersectParams::new(
        Intersection::ConstTheta,
        1.0,
        10000,
        Directionality::default(),
    );
    particle.intersect(machine, &intersect_params, &solver_params);
    assert!(matches!(
        particle.integration_status(),
        IntegrationStatus::IntersectedTimedOut
    ));
}
