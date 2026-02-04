//! Test NcGeometry functionality.

use std::path::PathBuf;

use dexter_equilibrium::{Geometry, NcFluxState, NcGeometryBuilder};
use ndarray::{Array1, Array2};
use rsl_interpolation::{Accelerator, Cache};

#[test]
#[allow(unused_variables)]
fn nc_geometry() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ1d = "steffen";
    let typ2d = "bicubic";
    let builder = NcGeometryBuilder::new(&path, typ1d, typ2d);
    let geometry = dbg!(builder.build().unwrap());

    let path: PathBuf = geometry.path();
    let typ1d: String = geometry.typ1d();
    let typ2d: String = geometry.typ2d();
    let baxis: f64 = geometry.baxis();
    let raxis: f64 = geometry.raxis();
    let zaxis: f64 = geometry.zaxis();
    let rgeo: f64 = geometry.rgeo();
    let r_wall: f64 = geometry.r_wall().unwrap();
    let shape: (usize, usize) = geometry.shape();
    let psi_state: NcFluxState = geometry.psi_state();
    let psip_state: NcFluxState = geometry.psip_state();
    let psi_wall: f64 = geometry.psi_wall().unwrap();
    let psip_wall: f64 = geometry.psip_wall().unwrap();
    let psi_array: Array1<f64> = geometry.psi_array();
    let psip_array: Array1<f64> = geometry.psip_array();
    let theta_array: Array1<f64> = geometry.theta_array();
    let r_array: Array1<f64> = geometry.r_array();
    let rlab_array: Array2<f64> = geometry.rlab_array();
    let zlab_array: Array2<f64> = geometry.rlab_array();
    let jacobian_array: Array2<f64> = geometry.rlab_array();

    let mut r_acc = Accelerator::new();
    let mut psi_acc = Accelerator::new();
    let mut psip_acc = Accelerator::new();
    let mut theta_acc = Accelerator::new();
    let mut cache = Cache::<f64>::new();
    let r = 0.2;
    let psi = 0.01;
    let psip = 0.015;
    let theta = 3.14;

    let psip_of_psi = geometry.psip_of_psi(psi, &mut psi_acc).unwrap();
    let psi_of_psip = geometry.psi_of_psip(psip, &mut psip_acc).unwrap();
    let r_of_psi = geometry.r_of_psi(psi, &mut psi_acc).unwrap();
    let r_of_psip = geometry.r_of_psip(psip, &mut psip_acc).unwrap();
    let psi_of_r = geometry.psi_of_r(r, &mut r_acc).unwrap();
    let psip_of_r = geometry.psip_of_r(r, &mut r_acc).unwrap();
    let rlab_of_psi = geometry
        .rlab_of_psi(psi, theta, &mut psi_acc, &mut theta_acc, &mut cache)
        .unwrap();
    let rlab_of_psip = geometry
        .rlab_of_psip(psip, theta, &mut psip_acc, &mut theta_acc, &mut cache)
        .unwrap();
    let zlab_of_psi = geometry
        .zlab_of_psi(psi, theta, &mut psi_acc, &mut theta_acc, &mut cache)
        .unwrap();
    let zlab_of_psip = geometry
        .zlab_of_psip(psip, theta, &mut psip_acc, &mut theta_acc, &mut cache)
        .unwrap();
    let jacobian_of_psi = geometry
        .jacobian_of_psi(psi, theta, &mut psi_acc, &mut theta_acc, &mut cache)
        .unwrap();
    let jacobian_of_psip = geometry
        .jacobian_of_psip(psip, theta, &mut psip_acc, &mut theta_acc, &mut cache)
        .unwrap();
}
