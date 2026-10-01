//! Calculation of the Energy in a 2D grid.

#![expect(clippy::min_ident_chars, reason = "Hamiltonian terms")]

use ndarray::{Array2, ArrayRef1, ArrayRef2, ArrayView2};
use rsl_interpolation::Accelerator2d;
use std::{f64::consts::TAU, mem::MaybeUninit};

use dexter_machine::{Machine, MagneticFlux::*};

/// Creates a `(θ, ψ/ψp)` meshgrid from the two 1D arrays.
///
/// This function is useful for creating the meshgrid needed for the [`energy_of_psi_grid`] and
/// [`energy_of_psip_grid`] functions.
///
/// To avoid extra allocations, this function returns two [`ArrayView2`]s.
///
/// # Indexing
///
/// The indexing order is [`XY`][ndarray::MeshIndex::XY], meaning:
/// - `theta_array` repeats on the rows of the grid.
/// - `flux_array` repeats on the columns of the grid.
///
/// Input:
/// ```text
///     `theta_array`:   | θ0, θ1, ..., θM|
///     `psi_array`:     | ψ0, ψ1, ..., ψN|
/// ```
///
/// Output:
///
/// ```text
///     `theta_grid`:    | θ0, θ1, ..., θM|
///                      | θ0, θ1, ..., θM|
///                      | θ0, θ1, ..., θM|
///                      | θ0, θ1, ..., θM|
///
///     `psi_grid`:      | ψ0, ψ0, ..., ψ0|
///                      | ψ1, ψ1, ..., ψ1|
///                      | ⋮,  ⋮,  ..., ⋮ |
///                      | ψN, ψN, ..., ψN|
/// ```
///
/// The energy grid would then be:
///
/// ```text
///     | E(θ0,ψ0), E(θ0,ψ1), ..., E(θ0,ψN) |
///     | E(θ1,ψ0), E(θ1,ψ1), ..., E(θ1,ψN) |
///
///     | E(θM,ψ0), E(θM,ψ1), ..., E(θM,ψN) |
/// ```
#[must_use]
pub fn create_poloidal_grid<'input_arrays>(
    theta_array: &'input_arrays ArrayRef1<f64>,
    flux_array: &'input_arrays ArrayRef1<f64>,
) -> (
    ArrayView2<'input_arrays, f64>,
    ArrayView2<'input_arrays, f64>,
) {
    ndarray::meshgrid((theta_array, flux_array), ndarray::MeshIndex::XY)
}

/// Calculates the energy on a 2D meshgrid of the `θ` and `ψ` arrays, in Normalized Units.
///
/// If any of the evaluations fail, the corresponding grid point is set to `f64::NAN`.
///
/// # Panics
///
/// This function panics if the two arrays do not have the same shape.
#[must_use]
pub fn energy_of_psi_grid(
    machine: Machine,
    pzeta: f64,
    mu: f64,
    theta_grid: &ArrayRef2<f64>,
    psi_grid: &ArrayRef2<f64>,
) -> Array2<f64> {
    assert_eq!(
        theta_grid.shape(),
        psi_grid.shape(),
        "Arrays must have the same shape"
    );

    let theta_array = theta_grid.row(0);
    let psi_array = psi_grid.column(0);
    let shape = (psi_array.len(), theta_array.len());
    let mut grid = Array2::<f64>::uninit(shape);

    assert_eq!(grid.nrows(), psi_array.len(), "sanity check");
    assert_eq!(grid.ncols(), theta_array.len(), "sanity check");

    let acc = &mut Accelerator2d::new();

    // Iterate though `psi_array` first to avoid unnecessarily recalculating `rho` and `g`.
    // Use the `Toroidal` flux for evaluation methods to avoid error propagation through double
    // interpolations.
    for m in 0..grid.nrows() {
        let psi = Toroidal(psi_array[m]);
        let Ok(psip) = machine.qfactor().eval_other(psi, acc.xacc()) else {
            std::hint::cold_path();
            nan_fill_row(&mut grid, m);
            continue;
        };
        let Ok(g) = machine.current().eval_g(psi, acc.xacc()) else {
            std::hint::cold_path();
            nan_fill_row(&mut grid, m);
            continue;
        };
        let rho = (pzeta + psip.value()) / g;
        for n in 0..grid.ncols() {
            let theta = theta_array[n].rem_euclid(TAU);
            if let Ok(b) = machine.bfield().eval_b(psi, theta, acc) {
                grid[[m, n]] = MaybeUninit::new((rho * b).powi(2) / 2.0 + mu * b)
            } else {
                std::hint::cold_path();
                grid[[m, n]] = MaybeUninit::new(f64::NAN)
            }
        }
    }

    // SAFETY: The loop passes from all elements and initializes them
    unsafe { grid.assume_init() }
}

/// Calculates the energy on a 2D meshgrid of the `θ` and `ψp` arrays, in Normalized Units.
///
/// If any of the evaluations fail, the corresponding grid point is set to `f64::NAN`.
///
/// # Panics
///
/// This function panics if the two arrays do not have the same shape.
#[must_use]
pub fn energy_of_psip_grid(
    machine: Machine,
    pzeta: f64,
    mu: f64,
    theta_grid: &ArrayRef2<f64>,
    psip_grid: &ArrayRef2<f64>,
) -> Array2<f64> {
    assert_eq!(
        theta_grid.shape(),
        psip_grid.shape(),
        "Arrays must have the same shape"
    );

    let theta_array = theta_grid.row(0);
    let psip_array = psip_grid.column(0);
    let shape = (psip_array.len(), theta_array.len());
    let mut grid = Array2::<f64>::uninit(shape);

    assert_eq!(grid.nrows(), psip_array.len(), "sanity check");
    assert_eq!(grid.ncols(), theta_array.len(), "sanity check");

    let acc = &mut Accelerator2d::new();

    // Iterate though `psip_array` first to avoid unnecessarily recalculating `rho`.
    for m in 0..grid.nrows() {
        let psip = Poloidal(psip_array[m]);
        let Ok(g) = machine.current().eval_g(psip, acc.xacc()) else {
            std::hint::cold_path();
            nan_fill_row(&mut grid, m);
            continue;
        };
        let rho = (pzeta + psip.value()) / g;
        for n in 0..grid.ncols() {
            let theta = theta_array[n].rem_euclid(TAU);
            if let Ok(b) = machine.bfield().eval_b(psip, theta, acc) {
                grid[[m, n]] = MaybeUninit::new((rho * b).powi(2) / 2.0 + mu * b)
            } else {
                std::hint::cold_path();
                grid[[m, n]] = MaybeUninit::new(f64::NAN)
            }
        }
    }

    // SAFETY: The loop passes from all elements and initializes them
    unsafe { grid.assume_init() }
}

/// Fills a 2D array's row with `f64::NAN`.
fn nan_fill_row(grid: &mut Array2<MaybeUninit<f64>>, row: usize) {
    grid.row_mut(row).fill(MaybeUninit::new(f64::NAN));
}

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::extract::TEST_NETCDF_PATH;
    use dexter_machine::*;
    use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    use ndarray::{Array1, arr1, arr2, array};

    #[test]
    fn poloidal_grid() {
        let theta_array = array![1.0, 2.0, 3.0, 4.0];
        let psi_array = array![10.0, 20.0, 30.0, 40.0];

        let (theta_grid, psi_grid) = dbg!(create_poloidal_grid(&theta_array, &psi_array));

        let expected_theta_grid = array![
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
        ];
        dbg!(&expected_theta_grid - &theta_grid);

        let expected_psi_grid = array![
            [10.0, 10.0, 10.0, 10.0],
            [20.0, 20.0, 20.0, 20.0],
            [30.0, 30.0, 30.0, 30.0],
            [40.0, 40.0, 40.0, 40.0],
        ];
        dbg!(&expected_psi_grid - &psi_grid);

        assert!(theta_grid.eq(&expected_theta_grid));
        assert!(psi_grid.eq(&expected_psi_grid));
    }

    #[test]
    fn gcmotion_check() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let theta_array = arr1(&vec![-1.0, 1.0]);
        let psi_array = arr1(&vec![0.01, 0.02]);

        let (theta_grid, psi_grid) = create_poloidal_grid(&theta_array, &psi_array);
        let energy_grid = energy_of_psi_grid(machine, -0.03, 1e-4, &theta_grid, &psi_grid);

        let expected = arr2(&[
            [0.00026296256388989685, 0.00026296256388989685],
            [0.00012897176092872728, 0.00012897176092872728],
        ])
        .to_owned();

        assert!(energy_grid.relative_eq(&expected, 1e-15, 1e-20));
    }

    #[test]
    fn toroidal_poloidal_equivalence() {
        let path = std::path::PathBuf::from(TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
        let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
        let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let acc = &mut Accelerator::new();
        // Avoid spline edges
        let theta_array = Array1::linspace(1e-3, TAU - 1e-3, 5);
        let psi_array = Array1::linspace(1e-3, 1.0 - 1e-3, 15) * qfactor.psi_last().value();
        let psip_array = psi_array.mapv(|psi| {
            let psi = Toroidal(psi);
            machine.qfactor().eval_other(psi, acc).unwrap().value()
        });

        let (theta_grid, psi_grid) = create_poloidal_grid(&theta_array, &psi_array);
        let (_, psip_grid) = create_poloidal_grid(&theta_array, &psip_array);

        let energy_of_psi = energy_of_psi_grid(machine, -0.03, 1e-4, &theta_grid, &psi_grid);
        let energy_of_psip = energy_of_psip_grid(machine, -0.03, 1e-4, &theta_grid, &psip_grid);

        dbg!(&energy_of_psi - &energy_of_psip);

        assert!(energy_of_psi.relative_eq(&energy_of_psip, 1e-7, 1e-7));
    }
}
