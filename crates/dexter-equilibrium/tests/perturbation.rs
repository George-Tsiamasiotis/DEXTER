use std::path::PathBuf;

use dexter_equilibrium::extract::TEST_NETCDF_PATH;
use dexter_equilibrium::{
    CosHarmonic, CosHarmonicCache, NcHarmonic, NcHarmonicBuilder, NcHarmonicCache, Perturbation,
};

#[test]
#[allow(unused_variables)]
#[rustfmt::skip]
fn empty_perturbation() {
    let p = Perturbation::zero();
    assert!(p.harmonics().is_empty());

    let caches: &mut Vec<CosHarmonicCache> = &mut p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    assert_eq!(p.p_of_psi(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.p_of_psip(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsi(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_dpsip(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dtheta(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dtheta(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dzeta(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dzeta(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psi_dt(psi, theta, zeta, t, caches).unwrap(), 0.0);
    assert_eq!(p.dp_of_psip_dt(psi, theta, zeta, t, caches).unwrap(), 0.0);
}

#[test]
#[allow(unused_variables)]
fn cos_perturbation() {
    let p = Perturbation::new(&[
        CosHarmonic::new(1e-3, 1, 1, 0.0),
        CosHarmonic::new(1e-3, 1, 2, 0.0),
        CosHarmonic::new(1e-3, 1, 3, 0.0),
    ]);
    let h0: CosHarmonic = p[0].clone();
    let h1: CosHarmonic = p[1].clone();
    let harmonics = p.harmonics();
    assert_eq!(p.count(), 3);

    let caches: &mut Vec<CosHarmonicCache> = &mut p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    p.p_of_psi(psi, theta, zeta, t, caches).unwrap();
    p.p_of_psip(psi, theta, zeta, t, caches).unwrap();
    p.dp_dpsi(psi, theta, zeta, t, caches).unwrap();
    p.dp_dpsip(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dtheta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dtheta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dzeta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dzeta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dt(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dt(psi, theta, zeta, t, caches).unwrap();
}

#[test]
#[allow(unused_variables)]
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

    let caches: &mut Vec<NcHarmonicCache> = &mut p.generate_caches();

    let (psi, theta, zeta, t) = (0.01, 1.0, 2.0, 0.0);
    p.p_of_psi(psi, theta, zeta, t, caches).unwrap();
    p.p_of_psip(psi, theta, zeta, t, caches).unwrap();
    p.dp_dpsi(psi, theta, zeta, t, caches).unwrap();
    p.dp_dpsip(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dtheta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dtheta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dzeta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dzeta(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psi_dt(psi, theta, zeta, t, caches).unwrap();
    p.dp_of_psip_dt(psi, theta, zeta, t, caches).unwrap();
}
