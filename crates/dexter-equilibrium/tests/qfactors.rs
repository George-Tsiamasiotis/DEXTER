//! Test Qfactors functionality.

use std::path::PathBuf;

use dexter_equilibrium::{
    EquilibriumType, FluxCommute, NcFluxState, NcQfactorBuilder, ParabolicQfactor, Qfactor,
    UnityQfactor,
};
use ndarray::Array1;
use rsl_interpolation::Accelerator;

#[test]
#[allow(unused_variables)]
fn nc_qfactor() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = NcQfactorBuilder::new(&path, typ);
    let qfactor = dbg!(builder.build().unwrap());

    let equilibrium_type: EquilibriumType = qfactor.equilibrium_type();
    let netcdf_version: semver::Version = qfactor.netcdf_version();
    let path: PathBuf = qfactor.path();
    let interp_type: String = qfactor.interp_type();
    let psi_state: NcFluxState = qfactor.psi_state();
    let psip_state: NcFluxState = qfactor.psip_state();
    let qaxis: f64 = qfactor.qaxis();
    let qwall: f64 = qfactor.qwall();
    let psi_wall: f64 = qfactor.psi_wall().unwrap();
    let psip_wall: f64 = qfactor.psip_wall().unwrap();
    let psi_array: Array1<f64> = qfactor.psi_array().unwrap();
    let psip_array: Array1<f64> = qfactor.psip_array().unwrap();
    let q_array: Array1<f64> = qfactor.q_array();

    let mut acc = Accelerator::new();
    let psi = 0.01;
    let psip = 0.015;
    let q_of_psi = qfactor.q_of_psi(psi, &mut acc).unwrap();
    let q_of_psip = qfactor.q_of_psip(psip, &mut acc).unwrap();
    let iota_of_psi = qfactor.iota_of_psi(psi, &mut acc).unwrap();
    let iota_of_psip = qfactor.iota_of_psip(psip, &mut acc).unwrap();
    let psip_of_psi = qfactor.psip_of_psi(psi, &mut acc).unwrap();
    let psi_of_psip = qfactor.psi_of_psip(psip, &mut acc).unwrap();
}

#[test]
#[allow(unused_variables)]
fn unity_qfactor() {
    let qfactor = dbg!(UnityQfactor::new());

    let mut acc = Accelerator::new();
    let p = 0.01;
    assert_eq!(qfactor.q_of_psi(p, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.q_of_psip(p, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.psip_of_psi(p, &mut acc).unwrap(), p);
    assert_eq!(qfactor.psi_of_psip(p, &mut acc).unwrap(), p);
    assert_eq!(qfactor.dpsip_dpsi(p, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.dpsi_dpsip(p, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.iota_of_psi(p, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.iota_of_psip(p, &mut acc).unwrap(), 1.0);
}

#[test]
#[allow(unused_variables)]
fn parabolic_qfactor() {
    let qfactor = dbg!(ParabolicQfactor::new(1.1, 3.8, 0.45));

    let qaxis: f64 = qfactor.qaxis();
    let qwall: f64 = qfactor.qwall();
    let psi_wall: f64 = qfactor.psi_wall();
    let psip_wall: f64 = qfactor.psip_wall();

    let mut acc = Accelerator::new();
    let psi = 0.01;
    let psip = 0.015;
    let q_of_psi = qfactor.q_of_psi(psi, &mut acc).unwrap();
    let q_of_psip = qfactor.q_of_psip(psip, &mut acc).unwrap();
    let iota_of_psi = qfactor.iota_of_psi(psi, &mut acc).unwrap();
    let iota_of_psip = qfactor.iota_of_psip(psip, &mut acc).unwrap();
    let psip_of_psi = qfactor.psip_of_psi(psi, &mut acc).unwrap();
    let psi_of_psip = qfactor.psi_of_psip(psip, &mut acc).unwrap();
}
