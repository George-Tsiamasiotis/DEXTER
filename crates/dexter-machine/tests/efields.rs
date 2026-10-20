//! Test Efields functionality.

#![allow(unused_variables)]

use dexter_machine::*;

#[test]
fn no_efield() {
    let efield = NoEfield;
    assert_eq!(efield.psi_state(), FluxCoordinateState::Good);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Good);

    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert_eq!(efield.eval_e(psi).unwrap(), 0.0);
    assert_eq!(efield.eval_e(psip).unwrap(), 0.0);
    assert_eq!(efield.eval_phi(psi).unwrap(), 0.0);
    assert_eq!(efield.eval_phi(psip).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv(psi).unwrap(), 0.0);
    assert_eq!(efield.eval_phi_deriv(psip).unwrap(), 0.0);
}

#[test]
fn toroidal_gaussian_efield() {
    let lcfs = MagneticFlux::Toroidal(0.45);

    let strength = 5e-4; // Normalized units
    let mean = lcfs * 0.9; // Peak close to the LCFS
    let width = 0.02 * lcfs.value(); // width of 2% of the ψ span.

    let efield = GaussianEfield::new(strength, mean, width);

    assert_eq!(efield.psi_state(), FluxCoordinateState::Good);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Bad);

    assert_eq!(efield.strength(), strength);
    assert_eq!(efield.mean(), mean);
    assert_eq!(efield.width(), width);

    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert!(efield.eval_e(psi).unwrap().is_finite());
    assert!(efield.eval_phi(psi).unwrap().is_finite());
    assert!(efield.eval_phi_deriv(psi).unwrap().is_finite());
    assert!(efield.eval_e(psip).is_err());
    assert!(efield.eval_phi(psip).is_err());
    assert!(efield.eval_phi_deriv(psip).is_err());
}

#[test]
fn poloidal_gaussian_efield() {
    let lcfs = MagneticFlux::Poloidal(0.45);

    let strength = 5e-4; // Normalized units
    let mean = lcfs * 0.9; // Peak close to the LCFS
    let width = 0.02 * lcfs.value(); // width of 2% of the ψ span.

    let efield = GaussianEfield::new(strength, mean, width);

    assert_eq!(efield.psi_state(), FluxCoordinateState::Bad);
    assert_eq!(efield.psip_state(), FluxCoordinateState::Good);

    assert_eq!(efield.strength(), strength);
    assert_eq!(efield.mean(), mean);
    assert_eq!(efield.width(), width);

    let psi = MagneticFlux::Toroidal(0.01);
    let psip = MagneticFlux::Poloidal(0.015);
    assert!(efield.eval_e(psip).unwrap().is_finite());
    assert!(efield.eval_phi(psip).unwrap().is_finite());
    assert!(efield.eval_phi_deriv(psip).unwrap().is_finite());
    assert!(efield.eval_e(psi).is_err());
    assert!(efield.eval_phi(psi).is_err());
    assert!(efield.eval_phi_deriv(psi).is_err());
}
