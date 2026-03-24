//! Test Perturbation's functionality.

#![allow(unused_variables)]

use std::path::PathBuf;

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::*;

#[test]
#[rustfmt::skip]
fn empty_perturbation() {
    let p = Perturbation::zero();
    assert!(p.harmonics().is_empty());

    let mut caches: Vec<CosHarmonicCache> = p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    assert_eq!(p.p_of_psi(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.p_of_psip(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsi(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsip(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dtheta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dtheta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dzeta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dzeta(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dt(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dt(psi, theta, zeta, t, &mut caches).unwrap(), 0.0);
}

#[test]
#[rustfmt::skip]
fn cos_toroidal_lcfs_perturbation() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let per = Perturbation::new(&[
        CosHarmonic::new(1e-3, lcfs, 1, 1, 0.0),
        CosHarmonic::new(1e-3, lcfs, 1, 2, 0.0),
        CosHarmonic::new(1e-3, lcfs, 1, 3, 0.0),
    ]);
    let h0: CosHarmonic = per[0].clone();
    let h1: CosHarmonic = per[1].clone();
    let harmonics = per.harmonics();
    assert_eq!(per.count(), 3);

    let mut c: Vec<CosHarmonicCache> = per.generate_caches();


    let (p, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);

    let _: f64 = per.p_of_psi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_dpsi(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psi_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(per.p_of_psip(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_dpsip(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psip_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero
}

#[test]
#[rustfmt::skip]
fn cos_poloidal_lcfs_perturbation() {
    let lcfs = LastClosedFluxSurface::Poloidal(0.45);
    let per = Perturbation::new(&[
        CosHarmonic::new(1e-3, lcfs, 1, 1, 0.0),
        CosHarmonic::new(1e-3, lcfs, 1, 2, 0.0),
        CosHarmonic::new(1e-3, lcfs, 1, 3, 0.0),
    ]);
    let h0: CosHarmonic = per[0].clone();
    let h1: CosHarmonic = per[1].clone();
    let harmonics = per.harmonics();
    assert_eq!(per.count(), 3);

    let mut c: Vec<CosHarmonicCache> = per.generate_caches();


    let (p, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);

    let _: f64 = per.p_of_psip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_dpsip(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dtheta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dzeta(p, theta, zeta, t, &mut c).unwrap();
    let _: f64 = per.dp_of_psip_dt(p, theta, zeta, t, &mut c).unwrap();

    assert!(per.p_of_psi(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_dpsi(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dtheta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dzeta(p, theta, zeta, t, &mut c).is_err());
    assert!(per.dp_of_psi_dt(p, theta, zeta, t, &mut c).is_ok()); // returns zero
}

#[test]
#[rustfmt::skip]
fn nc_perturbation() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let h1 = NcHarmonicBuilder::new(&path, "steffen", 2, 1)
        .build()
        .unwrap();
    let h2 = NcHarmonicBuilder::new(&path, "steffen", 2, 2)
        .build()
        .unwrap();
    let h3 = NcHarmonicBuilder::new(&path, "steffen", 3, 2)
        .build()
        .unwrap();

    let p = dbg!(Perturbation::new(&[h1, h2, h3]));
    let h0: NcHarmonic = p[0].clone();
    let h1: NcHarmonic = p[1].clone();
    let harmonics = p.harmonics();
    assert_eq!(harmonics.len(), 3);

    let mut caches: Vec<NcHarmonicCache> = p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    let _: f64 = p.p_of_psi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.p_of_psip(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_dpsi(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_dpsip(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psi_dtheta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psip_dtheta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psi_dzeta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psip_dzeta(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psi_dt(psi, theta, zeta, t, &mut caches).unwrap();
    let _: f64 = p.dp_of_psip_dt(psi, theta, zeta, t, &mut caches).unwrap();
}
