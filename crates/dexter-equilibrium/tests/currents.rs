//! Test NcCurrent functionality.

use std::path::PathBuf;

use dexter_equilibrium::{Current, NcCurrentBuilder, NcFluxState};
use ndarray::Array1;
use rsl_interpolation::Accelerator;

#[test]
#[allow(unused_variables)]
fn nc_current() {
    let path = PathBuf::from(dexter_equilibrium::extract::TEST_NETCDF_PATH);
    let typ = "steffen";
    let builder = NcCurrentBuilder::new(&path, typ);
    let current = dbg!(builder.build().unwrap());

    let path: PathBuf = current.path();
    let typ: String = current.typ();
    let len: usize = current.len();
    let psi_state: NcFluxState = current.psi_state();
    let psip_state: NcFluxState = current.psip_state();
    let psi_wall: f64 = current.psi_wall().unwrap();
    let psip_wall: f64 = current.psip_wall().unwrap();
    let psi_array: Array1<f64> = current.psi_array();
    let psip_array: Array1<f64> = current.psip_array();
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
