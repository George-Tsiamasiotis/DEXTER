//! Test the creation of an `EnergyPzetaPlane` from a set of `COMs`.

#![allow(unused_variables)]

use dexter_equilibrium::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
use dexter_equilibrium::*;
use dexter_simulate::*;
use parabola::Parabola;
use std::path::PathBuf;

#[test]
fn lar_energy_pzeta_parabola() {
    let lcfs = LastClosedFluxSurface::Toroidal(0.45);
    let qfactor = UnityQfactor::new(lcfs);
    let current = LarCurrent::new();
    let bfield = LarBfield::new();

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms
        .build_energy_pzeta_plane(&qfactor, &current, &bfield)
        .unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}

#[test]
fn toroidal_nc_energy_pzeta_parabola() {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms
        .build_energy_pzeta_plane(&qfactor, &current, &bfield)
        .unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}

#[test]
fn poloidal_nc_energy_pzeta_parabola() {
    let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
    let qfactor = NcQfactorBuilder::new(&path, "steffen").build().unwrap();
    let current = NcCurrentBuilder::new(&path, "steffen").build().unwrap();
    let bfield = NcBfieldBuilder::new(&path, "bicubic").build().unwrap();

    let coms = COMs {
        energy: None,
        pzeta: Some(-0.02),
        mu: Some(7e-6),
    };

    let plane = coms
        .build_energy_pzeta_plane(&qfactor, &current, &bfield)
        .unwrap();
    let axis: &Parabola = plane.axis_parabola();
    let left_wall: &Parabola = plane.left_wall_parabola();
    let right_wall: &Parabola = plane.right_wall_parabola();
    let tp_boundary: &TrappedPassingBoundary = plane.tp_boundary();

    assert_eq!(plane.mu(), coms.mu.unwrap());
}
