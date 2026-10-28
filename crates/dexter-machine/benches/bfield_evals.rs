//! Benchmark for the `Bfield` objects' evaluation methods.

#![allow(unused_results)]

use criterion::{Criterion, criterion_group, criterion_main};
use dexter_machine::extract::TEST_NETCDF_PATH;
use dexter_machine::*;

use std::hint::black_box;
use std::path::PathBuf;

fn bfield_evals(c: &mut Criterion) {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let acc = &mut Accelerator2d::new();
    let (psit, theta) = (MagneticFlux::Toroidal(0.1), 1.0);

    let lar_bfield = LarBfield;
    let nc_bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bicubic)
        .build()
        .unwrap();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield B(ψ, θ) evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_b(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_b(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield dB(ψ, θ)/dψ evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_deriv_flux(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_deriv_flux(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield dB(ψ, θ)/dθ evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_deriv_theta(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_deriv_theta(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield d²B(ψ, θ)/dψ² evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_deriv_flux2(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_deriv_flux2(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield d²B(ψ, θ)/dθ² evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_deriv_theta2(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_deriv_theta2(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Bfield d²B(ψ, θ)/dψdθ evaluation");

    group.bench_with_input("LarBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            lar_bfield
                .eval_deriv_mixed(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcBfield", &(psit, theta), |b, &(psit, theta)| {
        b.iter(|| {
            nc_bfield
                .eval_deriv_mixed(black_box(psit), black_box(theta), black_box(acc))
                .unwrap()
        });
    });
    group.finish();
}

criterion_group!(benches, bfield_evals);
criterion_main!(benches);
