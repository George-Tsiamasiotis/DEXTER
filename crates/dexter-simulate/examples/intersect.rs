//! `Particle::intersect` routine.

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;
use dexter_simulate::*;
use std::path::Path;

fn main() {
    analytical_equilibrium_intersect();
    numerical_equilibrium_intersect();
}

fn analytical_equilibrium_intersect() {
    // Equilibrium setup
    let qfactor = ParabolicQfactor::new(1.1, 3.8, LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let perturbation = Perturbation::new(&[
        CosHarmonic::new(1e-3, 2, 1, 0.0),
        CosHarmonic::new(1e-3, 3, 1, 0.0),
    ]);

    // Particle setup
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-4, 1e-6);
    let intersect_params = IntersectParams::new(Intersection::ConstTheta, 0.0, 100);
    let mut particle = Particle::new(&initial);

    // Calculate intersections
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );
    dbg!(&particle);
}

fn numerical_equilibrium_intersect() {
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
    let initial = InitialConditions::boozer(0.0, InitialFlux::Toroidal(0.2), 0.0, 0.0, 1e-6, 0.0);
    let intersect_params = IntersectParams::new(Intersection::ConstZeta, 0.0, 100);
    let mut particle = Particle::new(&initial);

    // Calculate intersections
    particle.intersect(
        &qfactor,
        &current,
        &bfield,
        &perturbation,
        &intersect_params,
        &SolverParams::default(),
    );
    dbg!(&particle);
}
