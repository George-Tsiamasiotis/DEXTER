//! Tests the creation of the [`StationaryCurve`] in multiple configurations.
//!
//! WARN: the order of the segments might change
//! TODO: find a way to make the segments reach the wall exactly

use std::f64::consts::PI;
use std::path::PathBuf;

use dexter_comspace::constants::SC_CONTOUR_FLUX_POINTS;
use dexter_comspace::*;
use dexter_machine::extract::{POLOIDAL_TEST_NETCDF_PATH, TOROIDAL_TEST_NETCDF_PATH};
use dexter_machine::*;

#[test]
fn lar_stationary_curve() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    check_lar_stationary_curve(curve, qfactor.psi_last().value());
}

#[test]
fn toroidal_lar_netcdf_stationary_curve() {
    let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = NcBfieldBuilder::new(
        &PathBuf::from(TOROIDAL_TEST_NETCDF_PATH),
        Interpolation2dType::Bicubic,
    )
    .build()
    .unwrap();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    check_lar_stationary_curve(curve, machine.qfactor().psi_last().value());
}

#[test]
fn poloidal_lar_netcdf_stationary_curve() {
    let qfactor = UnityQfactor::new(MagneticFlux::Poloidal(0.45));
    let current = LarCurrent::new();
    let bfield = NcBfieldBuilder::new(
        &PathBuf::from(POLOIDAL_TEST_NETCDF_PATH),
        Interpolation2dType::Bicubic,
    )
    .build()
    .unwrap();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    check_lar_stationary_curve(curve, machine.qfactor().psip_last().value());
}

fn check_lar_stationary_curve(curve: StationaryCurve, flux_last: f64) {
    assert_eq!(curve.segments.len(), 2);

    let psi_step = flux_last / SC_CONTOUR_FLUX_POINTS as f64;

    let seg1 = &curve.segments[0];
    let theta = seg1.theta();
    let flux = seg1.flux();

    assert!(theta.iter().all(|t| (*t - PI).abs() <= 1e-5));

    assert!(flux[0] >= flux_last - 1.5 * psi_step);
    assert!(flux[flux.len() - 1] < flux_last - 1.5 * psi_step);
    assert!((-&flux).iter().is_sorted());

    // ==========================

    let seg2 = curve.segments[1].clone();
    let theta = seg2.theta();
    let flux = seg2.flux();

    assert!(theta.iter().all(|t| t.abs() <= 1e-5));

    assert!(flux[0] < flux_last - 1.5 * psi_step);
    assert!(flux[flux.len() - 1] > flux_last - 1.5 * psi_step);
    assert!(flux.iter().is_sorted());
}
