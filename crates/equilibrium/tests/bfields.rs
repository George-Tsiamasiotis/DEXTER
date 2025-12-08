use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use rsl_interpolation::Accelerator;
use rsl_interpolation::Cache;
use std::path::PathBuf;

use equilibrium::Bfield;
use equilibrium::bfields;

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

    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let theta = 3.14;
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let mut cache = Cache::new();
    bfield
        .b(psip, theta, &mut xacc, &mut yacc, &mut cache)
        .unwrap();
    bfield
        .db_dpsip(psip, theta, &mut xacc, &mut yacc, &mut cache)
        .unwrap();
    bfield
        .db_dtheta(psip, theta, &mut xacc, &mut yacc, &mut cache)
        .unwrap();

    // =========================================

    assert_eq!(typ, bfield_typ);
    assert!(bfield_path.is_absolute());

    assert_eq!(psip_data.len(), bfield_shape.0);
    assert_eq!(theta_data.len(), bfield_shape.1);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(theta_data.ndim(), 1);
    assert_eq!(b_data.ndim(), 2);
}
