//! `Particle::integrate` routine.

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_simulate::*;
use std::path::Path;

fn main() {
    analytical_equilibrium_integration();
    numerical_equilibrium_integration();
}

fn analytical_equilibrium_integration() {
    // Equilibrium setup
    let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 4, 0.0),
    ]);

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 7e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e5),
        &SolverParams::default(),
    );
    dbg!(&particle);
}

fn numerical_equilibrium_integration() {
    // Equilibrium setup
    let path = Path::new("crates/dexter-simulate").join(TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();
    let perturbation = Perturbation::new(&[
        NcHarmonicBuilder::new(&path, "steffen", 2, 1)
            .with_phase_method(PhaseMethod::Interpolation)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .with_phase_method(PhaseMethod::Interpolation)
            .build()
            .unwrap(),
    ]);

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Poloidal(0.2), 0.0, 0.0, 1e-4, 7e-6);
    let mut particle = Particle::new(&initial);

    // Integrate
    particle.integrate(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        (0.0, 1e5),
        &SolverParams::default(),
    );
    dbg!(&particle);
}
