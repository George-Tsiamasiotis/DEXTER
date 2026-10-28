//! Test Perturbation's functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_machine::extract::TEST_NETCDF_PATH;
use dexter_machine::*;

#[test]
#[rustfmt::skip]
fn empty_perturbation() {
    let p = Perturbation::zero();

    let mut caches: DynModeCaches = p.generate_caches();

    let (psit, theta, zeta, t) = (MagneticFlux::Toroidal(0.01), 1.0, 2.0, 0.0);
    assert_eq!(p.eval_p(psit, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.eval_deriv_flux(psit, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.eval_deriv_theta(psit, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.eval_deriv_zeta(psit, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.eval_deriv_t(psit, theta, zeta, t, &mut caches).unwrap(), 0.0);
}

#[test]
#[rustfmt::skip]
fn cos_toroidal_lcfs_perturbation() {
    let lcfs = MagneticFlux::Toroidal(0.45);
    let per = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    ]);
    assert_eq!(per.count(), 3);

    let mut c: DynModeCaches = per.generate_caches();

    let (psit, theta, zeta, t) = (MagneticFlux::Toroidal(0.01), 1.0, 2.0, 0.0);

    assert!(per.eval_p(psit, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_flux(psit, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_theta(psit, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_zeta(psit, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_t(psit, theta, zeta, t, &mut c).unwrap().is_finite());

    let psip = MagneticFlux::Poloidal(0.015);
    assert!(per.eval_p(psip, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_flux(psip, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_theta(psip, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_zeta(psip, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_t(psip, theta, zeta, t, &mut c).is_ok()); // returns zero
}

#[test]
#[rustfmt::skip]
fn cos_poloidal_lcfs_perturbation() {
    let lcfs = MagneticFlux::Poloidal(0.45);
    let per = Perturbation::new(vec![
        Box::new(FluteMode::new(1e-3, lcfs, 1, 1, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
        Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    ]);
    assert_eq!(per.count(), 3);

    let mut c: DynModeCaches = per.generate_caches();


    let (psip, theta, zeta, t) = (MagneticFlux::Poloidal(0.01), 1.0, 2.0, 0.0);

    assert!(per.eval_p(psip, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_flux(psip, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_theta(psip, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_zeta(psip, theta, zeta, t, &mut c).unwrap().is_finite());
    assert!(per.eval_deriv_t(psip, theta, zeta, t, &mut c).unwrap().is_finite());

    let psit = MagneticFlux::Toroidal(0.01);
    assert!(per.eval_p(psit, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_flux(psit, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_theta(psit, theta, zeta, t, &mut c).is_err());
    assert!(per.eval_deriv_zeta(psit, theta, zeta, t, &mut c).is_err());
}

#[test]
#[rustfmt::skip]
fn nc_perturbation() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let typ = Interpolation1dType::Steffen;
    let m1 = NcFluteModeBuilder::new(&path, typ, 2, 1).build().unwrap();
    let m2 = NcFluteModeBuilder::new(&path, typ, 2, 2).build().unwrap();
    let m3 = NcFluteModeBuilder::new(&path, typ, 3, 2).build().unwrap();

    let per = dbg!(Perturbation::new(vec![Box::new(m1), Box::new(m2), Box::new(m3)]));
    assert_eq!(per.count(), 3);

    let mut caches: DynModeCaches = per.generate_caches();

    let (psit, theta, zeta, t) = (MagneticFlux::Toroidal(0.01), 1.0, 2.0, 0.0);
    assert!(per.eval_p(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_flux(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_theta(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_zeta(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_t(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
}

#[test]
#[rustfmt::skip]
fn mixed_perturbation() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let lcfs = MagneticFlux::Toroidal(0.45);
    let typ = Interpolation1dType::Cubic;
    let nc_mode = NcFluteModeBuilder::new(&path, typ, 2, 1).build().unwrap();
    let mode = FluteMode::new(1e-3, lcfs, 1, 3, 0.0);

    let per = dbg!(Perturbation::new(vec![Box::new(nc_mode), Box::new(mode)]));
    assert_eq!(per.count(), 2);

    let mut caches: DynModeCaches = per.generate_caches();

    let (psit, theta, zeta, t) = (MagneticFlux::Toroidal(0.01), 1.0, 2.0, 0.0);
    assert!(per.eval_p(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_flux(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_theta(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_zeta(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
    assert!(per.eval_deriv_t(psit, theta, zeta, t, &mut caches).unwrap().is_finite());
}
