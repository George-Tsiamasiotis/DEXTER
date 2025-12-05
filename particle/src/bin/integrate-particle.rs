use std::path::PathBuf;

use equilibrium::{Bfield, Currents, Harmonic, Perturbation, Qfactor};
use particle::*;

fn main() {
    let path = PathBuf::from("./data.nc");
    let qfactor = Qfactor::from_dataset(&path, "akima").unwrap();
    let currents = Currents::from_dataset(&path, "akima").unwrap();
    let bfield = Bfield::from_dataset(&path, "bicubic").unwrap();
    let harmonics = vec![
        Harmonic::from_dataset(&path, "akima", 1, 7).unwrap(),
        Harmonic::from_dataset(&path, "akima", 1, 9).unwrap(),
    ];
    let perturbation = Perturbation::from_harmonics(&harmonics);

    let initial = InitialConditions {
        time0: 0.0,
        theta0: 1.0,
        psip0: 0.002,
        rho0: 0.005,
        zeta0: 0.0,
        mu: 0.0,
    };

    let mut particle = Particle::new(&initial);
    particle.integrate(&qfactor, &bfield, &currents, &perturbation, (0.0, 100000.0));
    eprintln!("{:?}", particle);
}
