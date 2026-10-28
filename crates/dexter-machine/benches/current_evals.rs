//! Benchmark for the `Current` objects' evaluation methods.

#![allow(unused_results)]

use criterion::{Criterion, criterion_group, criterion_main};
use dexter_machine::extract::TEST_NETCDF_PATH;
use dexter_machine::*;

use std::hint::black_box;
use std::path::PathBuf;

fn current_evals(c: &mut Criterion) {
    let path = PathBuf::from(TEST_NETCDF_PATH);
    let acc = &mut Accelerator::new();
    let psit = MagneticFlux::Toroidal(0.01);

    let lar_current = LarCurrent;
    let nc_current = NcCurrentBuilder::new(&path, Interpolation1dType::Cubic)
        .build()
        .unwrap();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current g(ψ) evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| lar_current.eval_g(black_box(psit), black_box(acc)).unwrap());
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| nc_current.eval_g(black_box(psit), black_box(acc)).unwrap());
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current I(ψ) evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| lar_current.eval_i(black_box(psit), black_box(acc)).unwrap());
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| nc_current.eval_i(black_box(psit), black_box(acc)).unwrap());
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current dg(ψ)/dψ evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| {
            lar_current
                .eval_g_deriv(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| {
            nc_current
                .eval_g_deriv(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current dI(ψ)/dψ evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| {
            lar_current
                .eval_i_deriv(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| {
            nc_current
                .eval_i_deriv(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current d²g(ψ)/dψ² evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| {
            lar_current
                .eval_g_deriv2(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| {
            nc_current
                .eval_g_deriv2(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.finish();

    // ===========================================================================================

    let mut group = c.benchmark_group("Current d²I(ψ)/dψ² evaluation");

    group.bench_with_input("LarCurrent", &psit, |b, &psit| {
        b.iter(|| {
            lar_current
                .eval_i_deriv2(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.bench_with_input("NcCurrent", &psit, |b, &psit| {
        b.iter(|| {
            nc_current
                .eval_i_deriv2(black_box(psit), black_box(acc))
                .unwrap()
        });
    });
    group.finish();
}

criterion_group!(benches, current_evals);
criterion_main!(benches);
