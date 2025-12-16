use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use rsl_interpolation::Accelerator;
use std::path::PathBuf;

use equilibrium::*;

#[test]
fn test_nc_harmonic() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n);
    let harmonic = builder.build().unwrap();

    println!("{harmonic:?}");

    let harmonic_path = harmonic.path();
    let harmonic_typ = harmonic.typ();
    let harmonic_len = harmonic.len();
    let harmonic_m = harmonic.m();
    let harmonic_n = harmonic.n();
    let harmonic_phase_average = harmonic.phase_average();
    let harmonic_phase_resonance = harmonic.phase_resonance();
    let harmonic_phase_method = harmonic.phase_method();

    let psip_data = harmonic.psip_data();
    let a_data = harmonic.a_data();
    let phase_data = harmonic.phase_data();

    // =========================================

    assert!(harmonic_path.is_absolute());
    assert_eq!(typ, harmonic_typ);
    assert_eq!(m, harmonic_m);
    assert_eq!(n, harmonic_n);
    // The (1,2) resonance is inside the wall here
    assert!(harmonic_phase_average.is_none());
    assert!(harmonic_phase_resonance.is_some());
    assert!(matches!(harmonic_phase_method, PhaseMethod::Resonance));

    assert_eq!(psip_data.len(), harmonic_len);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(a_data.ndim(), 1);
    assert_eq!(phase_data.ndim(), 1);
}

#[test]
fn test_nc_harmonic_evals() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n);
    let harmonic = builder.build().unwrap();

    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let theta = 3.14;
    let zeta = 6.28;
    let mut acc = Accelerator::new();
    let mut cache = HarmonicCache::new();

    println!("{:?}", cache);

    // Normal
    assert!(
        harmonic
            .h(psip, theta, zeta, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dpsip(psip, theta, zeta, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dtheta(psip, theta, zeta, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dzeta(psip, theta, zeta, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dt(psip, theta, zeta, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(harmonic.a(psip, &mut acc).unwrap().is_finite());
    assert!(harmonic.da_dpsip(psip, &mut acc).unwrap().is_finite());
    assert!(harmonic.phase(psip, &mut acc).unwrap().is_finite());

    // Big θ and ζ
    assert!(
        harmonic
            .h(psip, 10000.0, 20000.0, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dpsip(psip, 10000.0, 20000.0, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dtheta(psip, 10000.0, 20000.0, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dzeta(psip, 10000.0, 20000.0, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        harmonic
            .dh_dt(psip, 10000.0, 20000.0, &mut acc, &mut cache)
            .unwrap()
            .is_finite()
    );

    // Out of bounds
    assert!(matches!(
        harmonic.a(100000.0, &mut acc),
        Err(EqError::DomainError(..))
    ));

    assert_eq!(cache.hits(), 6);
    assert_eq!(cache.misses(), 2);
}

#[test]
fn test_nc_harmonic_zero_phase_method() {
    use PhaseMethod::*;
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Zero);
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let mut acc = Accelerator::new();

    assert_eq!(harmonic.phase(0.1 * psip_wall, &mut acc).unwrap(), 0.0);
    assert_eq!(harmonic.phase(0.5 * psip_wall, &mut acc).unwrap(), 0.0);
    assert_eq!(harmonic.phase(0.9 * psip_wall, &mut acc).unwrap(), 0.0);
}

#[test]
fn test_nc_harmonic_average_phase_method() {
    use PhaseMethod::*;
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Average);
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let mut acc = Accelerator::new();

    let average = dbg!(harmonic.phase_average().unwrap());
    assert_eq!(harmonic.phase(0.1 * psip_wall, &mut acc).unwrap(), average);
    assert_eq!(harmonic.phase(0.5 * psip_wall, &mut acc).unwrap(), average);
    assert_eq!(harmonic.phase(0.9 * psip_wall, &mut acc).unwrap(), average);
}

#[test]
fn test_nc_harmonic_resonance_phase_method() {
    use PhaseMethod::*;
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Resonance);
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let mut acc = Accelerator::new();

    // WARN: Update this if `lar_netcdf.py` changes
    let must_be = -0.1790858645211414;
    dbg!(harmonic.phase_resonance().unwrap());
    dbg!(must_be);
    assert!((harmonic.phase(0.1 * psip_wall, &mut acc).unwrap() - must_be).abs() <= 1e-10);
    assert!((harmonic.phase(0.5 * psip_wall, &mut acc).unwrap() - must_be).abs() <= 1e-10);
    assert!((harmonic.phase(0.9 * psip_wall, &mut acc).unwrap() - must_be).abs() <= 1e-10);
}

#[test]
fn test_nc_harmonic_interpolation_phase_method() {
    use PhaseMethod::*;
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Interpolation);
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let mut acc = Accelerator::new();

    assert_ne!(
        harmonic.phase(0.1 * psip_wall, &mut acc).unwrap(),
        harmonic.phase(0.9 * psip_wall, &mut acc).unwrap()
    );
}

#[test]
fn test_nc_harmonic_custom_phase_method() {
    use PhaseMethod::*;
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let m = 2;
    let n = 1;
    let custom = 1.57;
    let builder = NcHarmonicBuilder::new(&path, typ, m, n).with_phase_method(Custom(custom));
    let harmonic = builder.build().unwrap();
    let psip_data = harmonic.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let mut acc = Accelerator::new();

    assert_eq!(harmonic.phase(0.1 * psip_wall, &mut acc).unwrap(), custom);
    assert_eq!(harmonic.phase(0.5 * psip_wall, &mut acc).unwrap(), custom);
    assert_eq!(harmonic.phase(0.9 * psip_wall, &mut acc).unwrap(), custom);
}
