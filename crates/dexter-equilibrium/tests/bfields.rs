//! Test Bfields functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::{
    Bfield, EquilibriumType, EvalError, LarBfield, NcBfieldBuilder, NcFluxState,
};
use ndarray::{Array1, Array2};
use rsl_interpolation::{Accelerator, Cache};

#[test]
#[rustfmt::skip]
fn lar_bfield() {
    let bfield = dbg!(LarBfield::new());

    let psi_acc = &mut Accelerator::new();
    let psip_acc = &mut Accelerator::new();
    let theta_acc = &mut Accelerator::new();
    let cache = &mut Cache::<f64>::new();
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = bfield.b_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_dpsi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_of_psi_dtheta(psi, theta, psi_acc, theta_acc, cache).unwrap();

    assert!(matches!(
        bfield.b_of_psip(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        bfield.db_dpsip(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        bfield.db_of_psip_dtheta(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
}

#[test]
#[rustfmt::skip]
fn nc_bfield() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let interp_type = "bicubic";
    let builder = NcBfieldBuilder::new(&path, interp_type);
    let bfield = dbg!(builder.build().unwrap());

    let equilibrium_type: EquilibriumType = bfield.equilibrium_type();
    let netcdf_version: semver::Version = bfield.netcdf_version();
    let path: PathBuf = bfield.path();
    let interp_type: String = bfield.interp_type();
    let baxis: f64 = bfield.baxis();
    let shape: (usize, usize) = bfield.shape();
    let psi_state: NcFluxState = bfield.psi_state();
    let psip_state: NcFluxState = bfield.psip_state();
    let psi_wall: f64 = bfield.psi_wall().unwrap();
    let psip_wall: f64 = bfield.psip_wall().unwrap();
    let psi_array: Array1<f64> = bfield.psi_array().unwrap();
    let psip_array: Array1<f64> = bfield.psip_array().unwrap();
    let theta_array: Array1<f64> = bfield.theta_array();
    let b_array: Array2<f64> = bfield.b_array();

    let psi_acc = &mut Accelerator::new();
    let psip_acc = &mut Accelerator::new();
    let theta_acc = &mut Accelerator::new();
    let cache = &mut Cache::<f64>::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = bfield.b_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.b_of_psip(psip, theta, psip_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_dpsi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_dpsip(psip, theta, psip_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_of_psi_dtheta(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = bfield.db_of_psip_dtheta(psip, theta, psip_acc, theta_acc, cache).unwrap();
}
