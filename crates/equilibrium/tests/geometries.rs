use equilibrium::extract::STUB_TEST_NETCDF_PATH;
use std::path::PathBuf;

use equilibrium::Geometry;
use equilibrium::geometries;

#[test]
fn test_nc_geometry() {
    let path = PathBuf::from(STUB_TEST_NETCDF_PATH);
    let typ1d = "cubic";
    let typ2d = "bicubic";
    let builder = geometries::NcGeometryBuilder::new(&path, typ1d, typ2d);
    let geometry = builder.build().unwrap();

    println!("{geometry:?}");

    let geometry_path = geometry.path();
    let geometry_typ1d = geometry.typ1d();
    let geometry_typ2d = geometry.typ2d();
    let geometry_shape = geometry.shape();
    let geometry_baxis = geometry.baxis();
    let geometry_raxis = geometry.raxis();
    let geometry_zaxis = geometry.zaxis();
    let geometry_rgeo = geometry.rgeo();
    let geometry_r_wall = geometry.r_wall();
    let geometry_psip_wall = geometry.psip_wall();
    let geometry_psi_wall = geometry.psi_wall();

    let psip_data = geometry.psip_data();
    let theta_data = geometry.theta_data();
    let psi_data = geometry.psi_data();
    let r_data = geometry.r_data();
    let rlab_data = geometry.rlab_data();
    let zlab_data = geometry.zlab_data();
    let jacobian_data = geometry.jacobian_data();

    let psip = 0.5 * geometry_psip_wall;
    let theta = 3.14;
    let r_from_psip = geometry.r(psip).unwrap();
    let psip_from_r = geometry.psip(r_from_psip).unwrap();
    geometry.rlab(psip, theta).unwrap();
    geometry.zlab(psip, theta).unwrap();

    // =========================================

    assert_eq!(typ1d, geometry_typ1d);
    assert_eq!(typ2d, geometry_typ2d);
    assert!(geometry_path.is_absolute());

    assert_eq!(psip_data.len(), geometry_shape.0);
    assert_eq!(theta_data.len(), geometry_shape.1);
    assert_eq!(psip_data.ndim(), 1);
    assert_eq!(theta_data.ndim(), 1);
    assert_eq!(psi_data.ndim(), 1);
    assert_eq!(r_data.ndim(), 1);
    assert_eq!(rlab_data.ndim(), 2);
    assert_eq!(zlab_data.ndim(), 2);
    assert_eq!(jacobian_data.ndim(), 2);

    assert_eq!(geometry_baxis, 1.5);
    assert_eq!(geometry_raxis, 1.75);
    assert_eq!(geometry_zaxis, 0.0);
    assert_eq!(geometry_rgeo, 1.7);

    assert!((psip_from_r - psip).abs() <= 1e-10);
    assert!((geometry_r_wall - geometry.r(geometry_psip_wall).unwrap()).abs() <= 1e-10);
    assert!((geometry_psi_wall - geometry.psi(geometry_psip_wall).unwrap()).abs() <= 1e-10);
}
