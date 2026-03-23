//! Test Geometry functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::{
    EquilibriumType, EvalError, FluxCommute, Geometry, LarGeometry, NcFluxState, NcGeometryBuilder,
};
use ndarray::{Array1, Array2};
use rsl_interpolation::{Accelerator, Cache};

#[test]
#[rustfmt::skip]
fn lar_geometry() {
    let geometry = dbg!(LarGeometry::new(2.0, 1.75, 0.5));

    let equilibrium_type: EquilibriumType = geometry.equilibrium_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let rlast: f64 = geometry.rlast();
    let psi_last: f64 = geometry.psi_last();
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let r_acc = &mut Accelerator::new();
    let psi_acc = &mut Accelerator::new();
    let psip_acc = &mut Accelerator::new();
    let theta_acc = &mut Accelerator::new();
    let cache = &mut Cache::<f64>::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = geometry.r_of_psi(psi, psi_acc).unwrap();
    let _: f64 = geometry.psi_of_r(r, r_acc).unwrap();
    let _: f64 = geometry.rlab_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.zlab_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();

    assert!(matches!(
        geometry.r_of_psip(psip, psip_acc),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.psip_of_r(r, r_acc),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.rlab_of_psip(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.zlab_of_psip(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.jacobian_of_psi(psi, theta, psi_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.jacobian_of_psip(psip, theta, psip_acc, theta_acc, cache),
        Err(EvalError::UndefinedEvaluation(..))
    ));
}

#[test]
#[rustfmt::skip]
fn nc_geometry() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ1d = "steffen";
    let typ2d = "bicubic";
    let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
    let geometry = dbg!(builder.build().unwrap());

    let equilibrium_type: EquilibriumType = geometry.equilibrium_type();
    let netcdf_version: semver::Version = geometry.netcdf_version();
    let path: PathBuf = geometry.path();
    let interp1d_type: String = geometry.interp1d_type();
    let interp2d_type: String = geometry.interp2d_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let zaxis: f64 = geometry.zaxis();
    let rgeo: f64 = geometry.rgeo();
    let rlast: f64 = geometry.rlast();
    let shape: (usize, usize) = geometry.shape();
    let psi_state: NcFluxState = geometry.psi_state();
    let psip_state: NcFluxState = geometry.psip_state();
    let psi_last: f64 = geometry.psi_last().unwrap();
    let psip_last: f64 = geometry.psip_last().unwrap();
    let psi_array: Array1<f64> = geometry.psi_array().unwrap();
    let psip_array: Array1<f64> = geometry.psip_array().unwrap();
    let theta_array: Array1<f64> = geometry.theta_array();
    let r_array: Array1<f64> = geometry.r_array();
    let rlab_array: Array2<f64> = geometry.rlab_array();
    let zlab_array: Array2<f64> = geometry.rlab_array();
    let jacobian_array: Array2<f64> = geometry.rlab_array();
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let r_acc = &mut Accelerator::new();
    let psi_acc = &mut Accelerator::new();
    let psip_acc = &mut Accelerator::new();
    let theta_acc = &mut Accelerator::new();
    let cache = &mut Cache::<f64>::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let _: f64 = geometry.psip_of_psi(psi, psi_acc).unwrap();
    let _: f64 = geometry.psi_of_psip(psip, psip_acc).unwrap();
    let _: f64 = geometry.r_of_psi(psi, psi_acc).unwrap();
    let _: f64 = geometry.r_of_psip(psip, psip_acc).unwrap();
    let _: f64 = geometry.psi_of_r(r, r_acc).unwrap();
    let _: f64 = geometry.psip_of_r(r, r_acc).unwrap();
    let _: f64 = geometry.rlab_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.rlab_of_psip(psip, theta, psip_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.zlab_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.zlab_of_psip(psip, theta, psip_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.jacobian_of_psi(psi, theta, psi_acc, theta_acc, cache).unwrap();
    let _: f64 = geometry.jacobian_of_psip(psip, theta, psip_acc, theta_acc, cache).unwrap();
}
