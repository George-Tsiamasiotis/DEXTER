//! Test Geometry functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_machine::*;
use ndarray::{Array1, ArrayView2};

#[test]
fn lar_geometry() {
    let geometry = dbg!(LarGeometry::new(2.0, 1.75, 0.5));

    assert_eq!(geometry.psit_state(), FluxCoordinateState::Good);
    assert_eq!(geometry.psip_state(), FluxCoordinateState::Bad);

    let machine_type: MachineType = geometry.machine_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let rlast: f64 = geometry.rlast();
    let psit_last: MagneticFlux = geometry.psit_last();
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let acc1 = &mut Accelerator::new();
    let acc2 = &mut Accelerator2d::new();
    let r = 0.2;
    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    let theta = 3.14;

    assert!(geometry.eval_r(psi, acc1).unwrap().is_finite());
    assert!(
        geometry
            .eval_psit_of_r(r, acc1)
            .unwrap()
            .value()
            .is_finite()
    );
    assert!(geometry.eval_rlab(psi, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_zlab(psi, theta, acc2).unwrap().is_finite());

    assert!(matches!(
        geometry.eval_r(psip, acc1),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.eval_psip_of_r(r, acc1),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.eval_rlab(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.eval_zlab(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.eval_jacobian(psi, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
    assert!(matches!(
        geometry.eval_jacobian(psip, theta, acc2),
        Err(EvalError::UndefinedEvaluation(..))
    ));
}

#[test]
#[rustfmt::skip]
fn nc_geometry() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ1d = Interpolation1dType::Cubic;
    let typ2d = Interpolation2dType::Bicubic;
    let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
    let geometry = dbg!(builder.build().unwrap());

    assert_eq!(geometry.psit_state(), FluxCoordinateState::Good);
    assert_eq!(geometry.psip_state(), FluxCoordinateState::Good);

    let machine_type: MachineType = geometry.machine_type();
    let netcdf_version: semver::Version = geometry.netcdf_version();
    let path: PathBuf = geometry.path();
    let interp1d_type: Interpolation1dType = geometry.interp1d_type();
    let interp2d_type: Interpolation2dType = geometry.interp2d_type();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let zaxis: f64 = geometry.zaxis();
    let rgeo: f64 = geometry.rgeo();
    let rlast: f64 = geometry.rlast();
    let shape: (usize, usize) = geometry.shape();
    let psit_state: FluxCoordinateState = geometry.psit_state();
    let psip_state: FluxCoordinateState = geometry.psip_state();
    let psit_last: MagneticFlux = geometry.psit_last().unwrap();
    let psip_last: MagneticFlux = geometry.psip_last().unwrap();
    let psit_array: &[f64] = geometry.psit_values().unwrap();
    let psip_array: &[f64] = geometry.psip_values().unwrap();
    let theta_array: &[f64] = geometry.theta_values();
    let r_array: &[f64] = geometry.r_values();
    let rlab_array: ArrayView2<f64> = geometry.rlab_array();
    let zlab_array: ArrayView2<f64> = geometry.rlab_array();
    let jacobian_array: ArrayView2<f64> = geometry.jacobian_array().unwrap();
    let rlab_last: Array1<f64> = geometry.rlab_last();
    let zlab_last: Array1<f64> = geometry.zlab_last();

    let acc1 = &mut Accelerator::new();
    let acc2 = &mut Accelerator2d::new();
    let r = 0.2;
    let psit = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    let theta = 3.14;

    assert!(geometry.eval_r(psit, acc1).unwrap().is_finite());
    assert!(geometry.eval_r(psip, acc1).unwrap().is_finite());
    assert!(geometry.eval_psit_of_r(r, acc1).unwrap().value().is_finite());
    assert!(geometry.eval_psip_of_r(r, acc1).unwrap().value().is_finite());
    assert!(geometry.eval_rlab(psit, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_rlab(psip, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_zlab(psit, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_zlab(psip, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_jacobian(psit, theta, acc2).unwrap().is_finite());
    assert!(geometry.eval_jacobian(psip, theta, acc2).unwrap().is_finite());
}
