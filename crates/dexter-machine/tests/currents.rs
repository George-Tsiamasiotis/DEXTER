//! Test Currents functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use approx::assert_abs_diff_eq;
use dexter_machine::*;

#[test]
fn lar_current() {
    let current = LarCurrent;
    assert_eq!(current.psit_state(), FluxCoordinateState::Good);
    assert_eq!(current.psip_state(), FluxCoordinateState::Good);

    let mut acc = Accelerator::new();
    let psit = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert_abs_diff_eq!(current.eval_g(psit, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(current.eval_g(psip, &mut acc).unwrap(), 1.0);
    assert_abs_diff_eq!(current.eval_i(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv2(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_g_deriv2(psip, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv2(psit, &mut acc).unwrap(), 0.0);
    assert_abs_diff_eq!(current.eval_i_deriv2(psip, &mut acc).unwrap(), 0.0);
}

#[test]
fn nc_current() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Cubic;
    let builder = NcCurrentBuilder::new(&path, typ);
    let current = dbg!(builder.build().unwrap());
    assert_eq!(current.psit_state(), FluxCoordinateState::Good);
    assert_eq!(current.psip_state(), FluxCoordinateState::Good);

    let machine_type: MachineType = current.machine_type();
    let netcdf_version: semver::Version = current.netcdf_version();
    let path: PathBuf = current.path();
    let interp_type: Interpolation1dType = current.interp_type();
    let psit_state: FluxCoordinateState = current.psit_state();
    let psip_state: FluxCoordinateState = current.psip_state();
    let psit_array: &[f64] = current.psit_values().unwrap();
    let psip_array: &[f64] = current.psip_values().unwrap();
    let g_array: &[f64] = current.g_values();
    let i_array: &[f64] = current.i_values();

    let acc = &mut Accelerator::new();
    let psit = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert!(current.eval_g(psit, acc).unwrap().is_finite());
    assert!(current.eval_g(psip, acc).unwrap().is_finite());
    assert!(current.eval_i(psit, acc).unwrap().is_finite());
    assert!(current.eval_i(psip, acc).unwrap().is_finite());
    assert!(current.eval_g_deriv(psit, acc).unwrap().is_finite());
    assert!(current.eval_g_deriv(psip, acc).unwrap().is_finite());
    assert!(current.eval_i_deriv(psit, acc).unwrap().is_finite());
    assert!(current.eval_i_deriv(psip, acc).unwrap().is_finite());
    assert!(current.eval_g_deriv2(psit, acc).unwrap().is_finite());
    assert!(current.eval_g_deriv2(psip, acc).unwrap().is_finite());
    assert!(current.eval_i_deriv2(psit, acc).unwrap().is_finite());
    assert!(current.eval_i_deriv2(psip, acc).unwrap().is_finite());
}
