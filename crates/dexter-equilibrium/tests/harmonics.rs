//! Test Harmonics' functionality.

#![allow(unused_variables)]

use dexter_equilibrium::*;
use ndarray::Array1;
use std::f64::consts::PI;
use std::path::PathBuf;

#[test]
#[rustfmt::skip]
fn cos_harmonic_toroidal_lcfs() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let har = dbg!(CosHarmonic::new(1e-3, lcfs, 3, 2, PI));
    assert_eq!(har.psi_state(), FluxCoordinateState::Good);
    assert_eq!(har.psip_state(), FluxCoordinateState::Bad);

    assert_eq!(har.equilibrium_type(), EquilibriumType::Analytical);
    assert_eq!(har.epsilon(), 1e-3);
    assert_eq!(har.m(), 3);
    assert_eq!(har.n(), 2);
    assert_eq!(har.phase(), PI);

    let p = 0.01;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = har.generate_cache();

    let _: f64 = har.alpha_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.phase_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.h_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_dpsi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psi_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psi_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psi_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(har.alpha_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(har.phase_of_psip(p, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(har.h_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_dpsip(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psip_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psip_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psip_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero

    assert_eq!(c.misses(), 1);
    assert_eq!(c.hits(), 4);
}

#[test]
#[rustfmt::skip]
fn cos_harmonic_poloidal_lcfs() {
    let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    let har = dbg!(CosHarmonic::new(1e-3, lcfs, 3, 2, PI));
    assert_eq!(har.psi_state(), FluxCoordinateState::Bad);
    assert_eq!(har.psip_state(), FluxCoordinateState::Good);

    assert_eq!(har.equilibrium_type(), EquilibriumType::Analytical);
    assert_eq!(har.epsilon(), 1e-3);
    assert_eq!(har.m(), 3);
    assert_eq!(har.n(), 2);
    assert_eq!(har.phase(), PI);

    let p = 0.01;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = har.generate_cache();

    let _: f64 = har.alpha_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.phase_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.h_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_dpsip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psip_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psip_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psip_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(har.alpha_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(har.phase_of_psi(p, theta, zeta, t, &mut c).is_ok()); // returns a constant
    assert!(har.h_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_dpsi(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psi_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psi_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(har.dh_of_psi_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero

    assert_eq!(c.misses(), 1);
    assert_eq!(c.hits(), 4);
}

#[test]
fn nc_harmonic() {
    let path = PathBuf::from(extract::TEST_NETCDF_PATH);
    let typ = "steffen";
    let (m, n) = (3, 2);
    use PhaseMethod::Interpolation;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Interpolation);
    let har = dbg!(builder.build().unwrap());
    assert_eq!(har.psi_state(), FluxCoordinateState::Good);
    assert_eq!(har.psip_state(), FluxCoordinateState::Good);

    assert_eq!(har.equilibrium_type(), EquilibriumType::Numerical);
    assert_eq!(har.m(), 3);
    assert_eq!(har.n(), 2);
    assert!(matches!(har.phase_method(), Interpolation));

    let netcdf_version: semver::Version = har.netcdf_version();
    let path: PathBuf = har.path();
    let interp_type: String = har.interp_type();
    let psi_state: FluxCoordinateState = har.psi_state();
    let psip_state: FluxCoordinateState = har.psip_state();
    let psi_last: f64 = har.psi_last().unwrap();
    let psip_last: f64 = har.psip_last().unwrap();
    let psi_array: Array1<f64> = har.psi_array().unwrap();
    let psip_array: Array1<f64> = har.psip_array().unwrap();
    let alpha_array: Array1<f64> = har.alpha_array();
    let phase_array: Array1<f64> = har.phase_array();

    let phase_average: Option<f64> = har.phase_average();
    let psi_phase_resonance: Option<f64> = har.psi_phase_resonance();
    let psip_phase_resonance: Option<f64> = har.psip_phase_resonance();

    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = har.generate_cache();

    let _: f64 = har.alpha_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.alpha_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.phase_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.phase_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.h_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.h_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psi_dtheta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psip_dtheta(psip, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psi_dzeta(psi, theta, zeta, t, &mut c).unwrap();
    let _: f64 = har.dh_of_psip_dzeta(psip, theta, zeta, t, &mut c).unwrap();
    assert_eq!(har.dh_of_psi_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(har.dh_of_psip_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
}
