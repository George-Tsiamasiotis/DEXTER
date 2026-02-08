//! Test Currents functionality.

use std::path::PathBuf;

use dexter_equilibrium::{Current, EquilibriumType, LarCurrent, NcCurrentBuilder, NcFluxState};
use ndarray::Array1;
use rsl_interpolation::Accelerator;

#[test]
#[allow(unused_variables)]
fn nc_current() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = NcCurrentBuilder::new(&path, typ);
    let current = dbg!(builder.build().unwrap());

    let equilibrium_type: EquilibriumType = current.equilibrium_type();
    let netcdf_version: semver::Version = current.netcdf_version();
    let path: PathBuf = current.path();
    let interp_type: String = current.interp_type();
    let psi_state: NcFluxState = current.psi_state();
    let psip_state: NcFluxState = current.psip_state();
    let psi_wall: f64 = current.psi_wall().unwrap();
    let psip_wall: f64 = current.psip_wall().unwrap();
    let psi_array: Array1<f64> = current.psi_array().unwrap();
    let psip_array: Array1<f64> = current.psip_array().unwrap();
    let g_array: Array1<f64> = current.g_array();
    let i_array: Array1<f64> = current.i_array();

    let mut acc = Accelerator::new();
    let psi = 0.01;
    let psip = 0.015;
    let g_of_psi = current.g_of_psi(psi, &mut acc).unwrap();
    let g_of_psip = current.g_of_psip(psip, &mut acc).unwrap();
    let i_of_psi = current.i_of_psi(psi, &mut acc).unwrap();
    let i_of_psip = current.i_of_psip(psip, &mut acc).unwrap();
    let dg_dpsi = current.dg_dpsi(psi, &mut acc).unwrap();
    let dg_dpsip = current.dg_dpsip(psip, &mut acc).unwrap();
    let di_dpsi = current.di_dpsi(psi, &mut acc).unwrap();
    let di_dpsip = current.di_dpsip(psip, &mut acc).unwrap();
}

#[test]
#[allow(unused_variables)]
fn lar_current() {
    let current = dbg!(LarCurrent::new());

    let mut acc = Accelerator::new();
    let p = 0.01;
    assert_eq!(current.g_of_psi(p, &mut acc).unwrap(), 1.0);
    assert_eq!(current.g_of_psip(p, &mut acc).unwrap(), 1.0);
    assert_eq!(current.i_of_psi(p, &mut acc).unwrap(), 0.0);
    assert_eq!(current.i_of_psip(p, &mut acc).unwrap(), 0.0);
    assert_eq!(current.dg_dpsi(p, &mut acc).unwrap(), 0.0);
    assert_eq!(current.dg_dpsip(p, &mut acc).unwrap(), 0.0);
    assert_eq!(current.di_dpsi(p, &mut acc).unwrap(), 0.0);
    assert_eq!(current.di_dpsip(p, &mut acc).unwrap(), 0.0);
}
