//! Test Efields functionality.

#![allow(unused_variables)]

use approx::assert_relative_eq;
use dexter_machine::*;

#[test]
fn no_efield() {
    let efield = NoEfield;
    assert_eq!(efield.psit_state(), FluxCoordinateState::Good);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Good);

    let psit = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert_eq!(efield.eval_e(psit).unwrap(), 0.0);
    assert_eq!(efield.eval_e(psip).unwrap(), 0.0);
    assert_eq!(efield.eval_phi(psit).unwrap(), 0.0);
    assert_eq!(efield.eval_phi(psip).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv(psit).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv(psip).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv2(psit).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv2(psip).unwrap(), 0.0);
}

#[test]
fn toroidal_gaussian_efield() {
    let lcfs = MagneticFlux::Toroidal(0.45);

    let strength = 5e-4; // Normalized units
    let mean = lcfs * 0.9; // Peak close to the LCFS
    let width = 0.02 * lcfs.value(); // width of 2% of the ψ span.

    let efield = GaussianEfield::new(strength, mean, width);

    assert_eq!(efield.psit_state(), FluxCoordinateState::Good);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Bad);

    assert_eq!(efield.strength(), strength);
    assert_eq!(efield.mean(), mean);
    assert_eq!(efield.width(), width);

    let psit = MagneticFlux::Toroidal(0.40);
    let psip = MagneticFlux::Poloidal(0.015);
    assert!(efield.eval_e(psit).unwrap().is_finite());
    assert!(efield.eval_phi(psit).unwrap().is_finite());
    assert!(efield.eval_phi_deriv(psit).unwrap().is_finite());
    assert!(efield.eval_phi_deriv2(psit).unwrap().is_finite());
    assert!(efield.eval_e(psip).is_err());
    assert!(efield.eval_phi(psip).is_err());
    assert!(efield.eval_phi_deriv(psip).is_err());
    assert!(efield.eval_phi_deriv2(psip).is_err());

    assert_relative_eq!(
        efield.eval_phi_deriv2(psit).unwrap(),
        (psit.value() - mean.value()) / width.powi(2) * efield.eval_e(psit).unwrap()
    );
}

#[test]
fn poloidal_gaussian_efield() {
    let lcfs = MagneticFlux::Poloidal(0.45);

    let strength = 5e-4; // Normalized units
    let mean = lcfs * 0.9; // Peak close to the LCFS
    let width = 0.02 * lcfs.value(); // width of 2% of the ψ span.

    let efield = GaussianEfield::new(strength, mean, width);

    assert_eq!(efield.psit_state(), FluxCoordinateState::Bad);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Good);

    assert_eq!(efield.strength(), strength);
    assert_eq!(efield.mean(), mean);
    assert_eq!(efield.width(), width);

    let psit = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert!(efield.eval_e(psip).unwrap().is_finite());
    assert!(efield.eval_phi(psip).unwrap().is_finite());
    assert!(efield.eval_phi_deriv(psip).unwrap().is_finite());
    assert!(efield.eval_phi_deriv2(psip).unwrap().is_finite());
    assert!(efield.eval_e(psit).is_err());
    assert!(efield.eval_phi(psit).is_err());
    assert!(efield.eval_phi_deriv(psit).is_err());
    assert!(efield.eval_phi_deriv2(psit).is_err());

    assert_relative_eq!(
        efield.eval_phi_deriv2(psip).unwrap(),
        (psip.value() - mean.value()) / width.powi(2) * efield.eval_e(psip).unwrap()
    );
}
