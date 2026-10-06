//! Tests the creation of the [`StationaryCurve`] in multiple configurations.
//!
//! WARN: the order of the segments might change
//! NOTE: segment starts/ends at the middle of the first/last two flux points

use approx::assert_relative_eq;
use dexter_machine::extract::POLOIDAL_TEST_NETCDF_PATH;
use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use dexter_comspace::constants::{SC_CONTOUR_FLUX_POINTS, SC_CONTOUR_THETA_POINTS};
use dexter_comspace::*;
use dexter_machine::{extract::TOROIDAL_TEST_NETCDF_PATH, *};

#[test]
fn lar_stationary_curve() {
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = LarBfield::new();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    assert_eq!(curve.segments.len(), 2);

    // ==========================

    let seg1 = curve.segments[0].clone();
    let theta = seg1.theta();
    let flux = seg1.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t - PI <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!((-&flux).iter().is_sorted());

    // ==========================

    let seg2 = curve.segments[1].clone();
    let theta = seg2.theta();
    let flux = seg2.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!(flux.iter().is_sorted());
}

#[test]
fn toroidal_lar_netcdf_stationary_curve() {
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.45));
    let current = LarCurrent::new();
    let bfield = NcBfieldBuilder::new(
        &PathBuf::from(TOROIDAL_TEST_NETCDF_PATH),
        Interpolation2dType::Bicubic,
    )
    .build()
    .unwrap();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    assert_eq!(curve.segments.len(), 2);

    // ==========================

    let seg1 = curve.segments[0].clone();
    let theta = seg1.theta();
    let flux = seg1.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t - PI <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!((-&flux).iter().is_sorted());

    // ==========================

    let seg2 = curve.segments[1].clone();
    let theta = seg2.theta();
    let flux = seg2.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psi_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!(flux.iter().is_sorted());
}

#[test]
fn poloidal_lar_netcdf_stationary_curve() {
    let qfactor = UnityQfactor::new(LastClosedFluxSurface::Poloidal(0.45));
    let current = LarCurrent::new();
    let bfield = NcBfieldBuilder::new(
        &PathBuf::from(POLOIDAL_TEST_NETCDF_PATH),
        Interpolation2dType::Bicubic,
    )
    .build()
    .unwrap();
    let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

    let curve = StationaryCurve::build(machine).unwrap();
    assert_eq!(curve.segments.len(), 2);

    // ==========================

    let seg1 = curve.segments[0].clone();
    let theta = seg1.theta();
    let flux = seg1.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t - PI <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psip_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psip_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!((-&flux).iter().is_sorted());

    // ==========================

    let seg2 = curve.segments[1].clone();
    let theta = seg2.theta();
    let flux = seg2.flux();

    assert!(
        theta
            .iter()
            .all(|t| *t <= TAU / SC_CONTOUR_THETA_POINTS as f64)
    );

    assert_relative_eq!(
        flux.last().copied().unwrap().abs(),
        machine.qfactor().psip_last().value() * (1.0 - (0.5 / SC_CONTOUR_FLUX_POINTS as f64))
    );
    assert_relative_eq!(
        flux.first().copied().unwrap().abs(),
        machine.qfactor().psip_last().value() / SC_CONTOUR_FLUX_POINTS as f64 / 2.0
    );
    assert!(flux.iter().is_sorted());
}
