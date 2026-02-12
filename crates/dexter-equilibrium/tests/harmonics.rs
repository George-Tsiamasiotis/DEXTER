//! Test Harmonics' functionality.

use dexter_equilibrium::{
    CosHarmonic, EquilibriumType, Harmonic, NcFluxState, NcHarmonicBuilder, PhaseMethod,
};
use ndarray::Array1;
use std::f64::consts::PI;
use std::path::PathBuf;

#[test]
#[allow(unused_variables)]
fn cos_harmonic() {
    let har = dbg!(CosHarmonic::new(1e-3, 3, 2, PI));

    assert_eq!(har.equilibrium_type(), EquilibriumType::Analytical);
    assert_eq!(har.ampl(), 1e-3);
    assert_eq!(har.m(), 3);
    assert_eq!(har.n(), 2);
    assert_eq!(har.phase(), PI);

    let psi = 0.01;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = har.get_default_cache();

    har.h_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    har.h_of_psip(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psi_dtheta(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psip_dtheta(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psi_dzeta(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psip_dzeta(psi, theta, zeta, t, &mut c).unwrap();

    assert_eq!(har.ampl_of_psi(psi, theta, zeta, t, &mut c).unwrap(), 1e-3);
    assert_eq!(har.ampl_of_psip(psi, theta, zeta, t, &mut c).unwrap(), 1e-3);
    assert_eq!(har.phase_of_psi(psi, theta, zeta, t, &mut c).unwrap(), PI);
    assert_eq!(har.phase_of_psip(psi, theta, zeta, t, &mut c).unwrap(), PI);
    assert_eq!(har.dh_dpsi(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(har.dh_dpsip(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(har.dh_of_psi_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(har.dh_of_psip_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);

    assert_eq!(c.hits(), 5); // not all evaluations use the cache
    assert_eq!(c.misses(), 1);
}

#[test]
#[allow(unused_variables)]
fn nc_harmonic() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ = "steffen";
    let (m, n) = (3, 2);
    use PhaseMethod::Interpolation;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Interpolation);
    let har = dbg!(builder.build().unwrap());

    assert_eq!(har.equilibrium_type(), EquilibriumType::Numerical);
    assert_eq!(har.m(), 3);
    assert_eq!(har.n(), 2);
    assert!(matches!(har.phase_method(), Interpolation));

    let netcdf_version: semver::Version = har.netcdf_version();
    let path: PathBuf = har.path();
    let interp_type: String = har.interp_type();
    let psi_state: NcFluxState = har.psi_state();
    let psip_state: NcFluxState = har.psip_state();
    let psi_wall: f64 = har.psi_wall().unwrap();
    let psip_wall: f64 = har.psip_wall().unwrap();
    let psi_array: Array1<f64> = har.psi_array().unwrap();
    let psip_array: Array1<f64> = har.psip_array().unwrap();
    let alpha_array: Array1<f64> = har.alpha_array();
    let phase_array: Array1<f64> = har.phase_array();

    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;
    let zeta = 1.0;
    let t = 8.0;
    let mut c = har.get_default_cache();

    har.ampl_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    har.ampl_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    har.phase_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    har.phase_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    har.h_of_psi(psi, theta, zeta, t, &mut c).unwrap();
    har.h_of_psip(psip, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psi_dtheta(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psip_dtheta(psip, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psi_dzeta(psi, theta, zeta, t, &mut c).unwrap();
    har.dh_of_psip_dzeta(psip, theta, zeta, t, &mut c).unwrap();
    assert_eq!(har.dh_of_psi_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
    assert_eq!(har.dh_of_psip_dt(psi, theta, zeta, t, &mut c).unwrap(), 0.0);
}
