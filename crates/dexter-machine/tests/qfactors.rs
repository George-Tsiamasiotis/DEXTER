//! Test Qfactors functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use approx::{assert_abs_diff_eq, assert_relative_eq};
use dexter_machine::{MagneticFlux::*, *};

#[test]
fn unity_qfactor() {
    let qfactor = dbg!(UnityQfactor::new(Toroidal(0.45)));

    assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let mut acc = Accelerator::new();
    let psit = Toroidal(0.01);
    let psip = Poloidal(0.015);
    assert_abs_diff_eq!(qfactor.eval_q(psit, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_q(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_q_deriv(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(qfactor.eval_q_deriv(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(qfactor.eval_iota(psit, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_iota(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_iota_deriv(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(qfactor.eval_iota_deriv(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(qfactor.eval_other(psit, &mut acc).unwrap(), Poloidal(0.01));
    assert_abs_diff_eq!(qfactor.eval_other(psip, &mut acc).unwrap(), Toroidal(0.015));
    assert_abs_diff_eq!(qfactor.eval_deriv_of_other(psit, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(qfactor.eval_deriv_of_other(psip, &mut acc).unwrap(), 1.0);

    assert!(qfactor.eval_psit_of_q(1.0, &mut acc).is_err());
    assert!(qfactor.eval_psip_of_q(1.0, &mut acc).is_err());

    assert_eq!(qfactor.psit_last(), Toroidal(0.45));
    assert_eq!(qfactor.psip_last(), Poloidal(0.45));
}

#[test]
fn parabolic_qfactor() {
    let qfactor = dbg!(ParabolicQfactor::new(1.1, 3.8, Toroidal(0.45)));

    assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let qaxis: f64 = qfactor.qaxis();
    let qlast: f64 = qfactor.qlast();
    let psit_last: MagneticFlux = qfactor.psit_last();
    let psip_last: MagneticFlux = qfactor.psip_last();

    let acc = &mut Accelerator::new();
    let q = 1.3;
    let psit = Toroidal(0.01);
    let psip = Poloidal(0.015);
    assert!(qfactor.eval_q(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_q(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_q_deriv(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_q_deriv(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota_deriv(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota_deriv(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_other(psit, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_other(psip, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_deriv_of_other(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_of_other(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_wrt_other(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_wrt_other(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_psit_of_q(q, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_psip_of_q(q, acc).unwrap().value().is_finite());
}

#[test]
fn nc_qfactor() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Cubic;
    let builder = NcQfactorBuilder::new(&path, typ);
    let qfactor = dbg!(builder.build().unwrap());

    assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
    assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

    let machine_type: MachineType = qfactor.machine_type();
    let netcdf_version: semver::Version = qfactor.netcdf_version();
    let path: PathBuf = qfactor.path();
    let interp_type: Interpolation1dType = qfactor.interp_type();
    let psit_state: FluxCoordinateState = qfactor.psit_state();
    let psip_state: FluxCoordinateState = qfactor.psip_state();
    let qaxis: f64 = qfactor.qaxis();
    let qlast: f64 = qfactor.qlast();
    let psi_last: MagneticFlux = qfactor.psit_last();
    let psip_last: MagneticFlux = qfactor.psip_last();
    let psi_array: &[f64] = qfactor.psit_values();
    let psip_array: &[f64] = qfactor.psip_values();
    let q_array: &[f64] = qfactor.q_values();

    let acc = &mut Accelerator::new();
    let q = 1.1;
    let psit = Toroidal(0.01);
    let psip = Poloidal(0.015);
    assert!(qfactor.eval_q(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_q(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_q_deriv(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_q_deriv(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota_deriv(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_iota_deriv(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_other(psit, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_other(psip, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_deriv_of_other(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_of_other(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_wrt_other(psit, acc).unwrap().is_finite());
    assert!(qfactor.eval_deriv_wrt_other(psip, acc).unwrap().is_finite());
    assert!(qfactor.eval_psit_of_q(q, acc).unwrap().value().is_finite());
    assert!(qfactor.eval_psip_of_q(q, acc).unwrap().value().is_finite());
}

#[test]
fn nc_qfactor_inverse_q() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Akima;
    let builder = NcQfactorBuilder::new(&path, typ);
    let qfactor = dbg!(builder.build().unwrap());
    let acc = &mut Accelerator::new();

    let psit = Toroidal(0.1);
    let q_of_psit = qfactor.eval_q(psit, acc).unwrap();
    let psit_inverse = qfactor.eval_psit_of_q(q_of_psit, acc).unwrap();
    assert_relative_eq!(psit, psit_inverse, epsilon = 1e-6);

    let psip = Poloidal(0.1);
    let q_of_psip = qfactor.eval_q(psip, acc).unwrap();
    let psip_inverse = qfactor.eval_psip_of_q(q_of_psip, acc).unwrap();
    assert_relative_eq!(psip, psip_inverse, epsilon = 1e-6);
}
