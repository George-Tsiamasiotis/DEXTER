use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use rsl_interpolation::Accelerator;
use std::path::PathBuf;

use equilibrium::Current;
use equilibrium::currents;

#[test]
fn test_nc_current() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = currents::NcCurrentBuilder::new(&path, typ);
    let current = builder.build().unwrap();

    println!("{current:?}");

    let current_path = current.path();
    let current_typ = current.typ();
    let current_len = current.len();

    let psip_data = current.psip_data();
    let g_data = current.g_data();
    let i_data = current.i_data();

    let psip_wall = psip_data.last().unwrap();
    let psip = 0.5 * psip_wall;
    let mut acc = Accelerator::new();
    current.g(psip, &mut acc).unwrap();
    current.i(psip, &mut acc).unwrap();
    current.dg_dpsip(psip, &mut acc).unwrap();
    current.di_dpsip(psip, &mut acc).unwrap();

    // =========================================

    assert_eq!(typ, current_typ);
    assert!(current_path.is_absolute());

    assert_eq!(psip_data.len(), current_len);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(g_data.ndim(), 1);
    assert_eq!(i_data.ndim(), 1);
}

#[test]
fn test_g_axis_value() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = currents::NcCurrentBuilder::new(&path, typ);
    let current = builder.build().unwrap();

    let mut acc = Accelerator::new();
    assert!(current.g(0.0, &mut acc).unwrap().abs() - 1.0 <= 1e-3);
}

#[test]
fn test_lar_current() {
    let current = currents::Lar;

    println!("{current:?}");

    let psip = 0.1;
    let mut acc = Accelerator::new();
    assert_eq!(current.g(psip, &mut acc).unwrap(), 1.0);
    assert_eq!(current.i(psip, &mut acc).unwrap(), 0.0);
    assert_eq!(current.dg_dpsip(psip, &mut acc).unwrap(), 0.0);
    assert_eq!(current.di_dpsip(psip, &mut acc).unwrap(), 0.0);
}
