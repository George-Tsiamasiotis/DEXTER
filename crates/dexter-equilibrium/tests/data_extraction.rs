//! Test netCDF extraction functions.

use dexter_equilibrium::extract::netcdf_fields::*;
use dexter_equilibrium::extract::*;
use ndarray::{Array1, Array2, Array3};

use std::path::PathBuf;

#[test]
#[allow(unused_variables)]
fn main() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let file = open(&path).unwrap();

    let baxis: f64 = extract_scalar(&file, BAXIS).unwrap();
    let psis: Array1<f64> = extract_1d_array(&file, PSI).unwrap();
    let b_norm: Array2<f64> = extract_2d_array(&file, B_NORM).unwrap();
    let _alphas: Array3<f64> = extract_3d_array(&file, ALPHAS).unwrap();
    let (alphas_norm_m2_n1, phases_m2_n1) = extract_harmonic_arrays::<f64>(&file, 2, 1).unwrap();

    let theta_var = extract_variable(&file, THETA).unwrap();
}
