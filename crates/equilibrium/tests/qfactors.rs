use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use ndarray::Array1;
use rsl_interpolation::Accelerator;
use std::path::PathBuf;

use equilibrium::Qfactor;
use equilibrium::qfactors;

#[test]
fn test_nc_qfactor() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = qfactors::NcQfactorBuilder::new(&path, typ);
    let qfactor = builder.build().unwrap();

    println!("{qfactor:?}");

    let qfactor_path = qfactor.path();
    let qfactor_typ = qfactor.typ();
    let qfactor_len = qfactor.len();

    let psip_data = qfactor.psip_data();
    let q_data = qfactor.q_data();
    let psi_data = qfactor.psi_data();

    // =========================================

    assert_eq!(typ, qfactor_typ);
    assert!(qfactor_path.is_absolute());

    assert_eq!(psip_data.len(), qfactor_len);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(q_data.ndim(), 1);
    assert_eq!(psi_data.ndim(), 1);
}

#[test]
fn test_nc_qfactor_evals() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = qfactors::NcQfactorBuilder::new(&path, typ);
    let qfactor = builder.build().unwrap();

    let psip_data = qfactor.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let mut acc = Accelerator::new();

    // Normal
    assert!(qfactor.q(psip, &mut acc).unwrap().is_finite());
    assert!(qfactor.psi(psip, &mut acc).unwrap().is_finite());
    assert!(qfactor.dpsi_dpsip(psip, &mut acc).unwrap().is_finite());

    // Out of bounds
    assert!(matches!(
        qfactor.q(100000.0, &mut acc),
        Err(equilibrium::EqError::DomainError(..))
    ));
}

#[test]
fn test_nc_qfactor_q_dpsidpsip() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = qfactors::NcQfactorBuilder::new(&path, typ);
    let qfactor = builder.build().unwrap();

    let mut acc = Accelerator::new();
    let psip_wall = qfactor.psip_data().last().copied().unwrap();
    let q_wall = qfactor.q(psip_wall, &mut acc).unwrap();
    // Do not go to close to the edges, since the interpolation might deviate a bit
    let psips = Array1::linspace(0.02 * psip_wall, 0.98 * psip_wall, 100);

    for psip in psips.iter() {
        assert!(
            (qfactor.q(*psip, &mut acc).unwrap() - qfactor.dpsi_dpsip(*psip, &mut acc).unwrap())
                .abs()
                < q_wall * 1e-4
        )
    }
}

#[test]
fn test_unity_qfactor() {
    let qfactor = qfactors::Unity;

    println!("{qfactor:?}");

    let psip = 0.1;
    let mut acc = Accelerator::new();
    assert_eq!(qfactor.q(psip, &mut acc).unwrap(), 1.0);
    assert_eq!(qfactor.psi(psip, &mut acc).unwrap(), psip);
    assert_eq!(qfactor.dpsi_dpsip(psip, &mut acc).unwrap(), 1.0);
}
