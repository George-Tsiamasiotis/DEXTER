use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use rsl_interpolation::Accelerator;
use std::path::PathBuf;

use equilibrium::*;

#[test]
fn test_nc_perturbation_no_harmonics() {
    let perturbation = NcPerturbation::from_harmonics(&[]);

    println!("{perturbation:?}");

    assert!(perturbation.is_empty());
    assert_eq!(perturbation.get_harmonics().len(), 0);
    assert_eq!(perturbation.len(), 0);

    let psip = 0.005;
    let theta = 3.14;
    let zeta = 6.28;
    let mut acc = Accelerator::new();
    let mut hcaches = vec![HarmonicCache::new()];

    assert_eq!(
        perturbation
            .p(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap(),
        0.0
    );
    assert_eq!(
        perturbation
            .dp_dpsip(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap(),
        0.0
    );
    assert_eq!(
        perturbation
            .dp_dtheta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap(),
        0.0
    );
    assert_eq!(
        perturbation
            .dp_dzeta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap(),
        0.0
    );
    assert_eq!(
        perturbation
            .dp_dt(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap(),
        0.0
    );
}

#[test]
fn test_nc_perturbation_one_harmonic() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let builder = NcHarmonicBuilder::new(&path, "steffen", 2, 1);
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let perturbation = NcPerturbation::from_harmonics(&[harmonic]);

    assert!(!perturbation.is_empty());
    println!("{perturbation:?}");

    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let theta = 3.14;
    let zeta = 6.28;
    let mut acc = Accelerator::new();
    let mut hcaches = vec![HarmonicCache::new(); perturbation.len()];

    // Normal
    assert!(
        perturbation
            .p(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dpsip(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dtheta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dzeta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dt(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );

    // Big θ and ζ
    assert!(
        perturbation
            .dp_dzeta(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .p(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dpsip(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dtheta(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dt(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );

    // Out of bounds
    assert!(matches!(
        dbg!(perturbation.dp_dtheta(10000.0, theta, zeta, &mut acc, &mut hcaches)),
        Err(equilibrium::EqError::DomainError(..))
    ));
}

#[test]
fn test_nc_perturbation_multiple_harmonics() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let perturbation = NcPerturbation::from_harmonics(&[
        NcHarmonicBuilder::new(&path, "steffen", 2, 1)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 2, 2)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 1)
            .build()
            .unwrap(),
        NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .build()
            .unwrap(),
    ]);
    let psip_data = perturbation.get_harmonics().first().unwrap().psip_data();

    assert!(!perturbation.is_empty());
    println!("{perturbation:?}");

    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let theta = 3.14;
    let zeta = 6.28;
    let mut acc = Accelerator::new();
    let mut hcaches = vec![HarmonicCache::new(); perturbation.len()];

    // Normal
    assert!(
        perturbation
            .dp_dtheta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .p(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dpsip(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dzeta(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dt(psip, theta, zeta, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );

    // Big θ and ζ
    assert!(
        perturbation
            .dp_dpsip(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .p(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dtheta(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dzeta(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );
    assert!(
        perturbation
            .dp_dt(psip, 10000.0, 20000.0, &mut acc, &mut hcaches)
            .unwrap()
            .is_finite()
    );

    // Out of bounds
    assert!(matches!(
        dbg!(perturbation.dp_dtheta(10000.0, theta, zeta, &mut acc, &mut hcaches)),
        Err(equilibrium::EqError::DomainError(..))
    ));
}
