use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use rsl_interpolation::Accelerator;
use rsl_interpolation::Cache;
use std::path::PathBuf;

use equilibrium::bfields;
use equilibrium::*;

#[test]
fn test_nc_bfield() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "bicubic";
    let builder = bfields::NcBfieldBuilder::new(&path, typ);
    let bfield = builder.build().unwrap();

    println!("{bfield:?}");

    let bfield_path = bfield.path();
    let bfield_typ = bfield.typ();
    let bfield_shape = bfield.shape();

    let psip_data = bfield.psip_data();
    let theta_data = bfield.theta_data();
    let b_data = bfield.b_data();

    let db_dpsip_data = bfield.db_dpsip_data();
    let db_dtheta_data = bfield.db_dtheta_data();

    // =========================================

    assert_eq!(typ, bfield_typ);
    assert!(bfield_path.is_absolute());

    assert_eq!(psip_data.len(), bfield_shape.0);
    assert_eq!(theta_data.len(), bfield_shape.1);
    assert_eq!(db_dpsip_data.shape(), &[bfield_shape.0, bfield_shape.1]);
    assert_eq!(db_dtheta_data.shape(), &[bfield_shape.0, bfield_shape.1]);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(theta_data.ndim(), 1);
    assert_eq!(b_data.ndim(), 2);
}

#[test]
fn test_nc_bfield_evals() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "bicubic";
    let builder = bfields::NcBfieldBuilder::new(&path, typ);
    let bfield = builder.build().unwrap();

    let psip_data = bfield.psip_data();
    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let theta = 3.14;
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let mut cache = Cache::new();

    // Normal
    assert!(
        bfield
            .b(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        bfield
            .db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        bfield
            .db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );

    // Big θ
    assert!(
        bfield
            .b(psip, 10000.0, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        bfield
            .db_dpsip(psip, 10000.0, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );
    assert!(
        bfield
            .db_dtheta(psip, 10000.0, &mut xacc, &mut yacc, &mut cache)
            .unwrap()
            .is_finite()
    );

    // Out of bounds
    assert!(matches!(
        bfield.b(100000.0, theta, &mut xacc, &mut yacc, &mut cache),
        Err(EqError::DomainError(..))
    ));
}
