//! Evaluation of the Hamiltonian in 1D arrays and 2D meshgrids.

#![expect(clippy::min_ident_chars, reason = "Hamiltonian terms")]

use ndarray::{Array1, Array2, ArrayRef1, ArrayRef2, ArrayView2};
use rsl_interpolation::Accelerator2d;
use std::{f64::consts::TAU, mem::MaybeUninit};

use dexter_machine::{FluxCoordinateState, Machine, MagneticFlux::*};

use crate::{COMError, GCHamiltonian};

impl GCHamiltonian {
    /// Evaluates the Hamiltonian in a `([θ], [ψ]; Pζ, μ)` 1D array.
    ///
    /// Use this method when you *know* that `psi` is a [`Good`] coordinate.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psi` is not a [`Good`] coordinate.
    /// + Returns a [`COMError::ShapeMismatch`] if `theta_array` and `psi_array` do not have the
    ///   same length.
    /// + Returns a [`COMError::EvalError`] if `psi_array` contains out-of-bounds values.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_comspace::*;
    /// # use ndarray::Array1;
    /// # use std::f64::consts::TAU;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 10);
    /// let psi_array = Array1::linspace(0.0, machine.psi_last().value(), 10);
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let energy_array =
    ///     GCHamiltonian::eval_psi_array(machine, &theta_array, &psi_array, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_psi_array(
        machine: Machine,
        theta_array: &ArrayRef1<f64>,
        psi_array: &ArrayRef1<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array1<f64>, COMError> {
        if machine.psi_state() != FluxCoordinateState::Good {
            return Err(COMError::BadFluxCoordinate);
        }
        if !(theta_array.len() == psi_array.len()) {
            return Err(COMError::ShapeMismatch);
        }

        let acc = &mut Accelerator2d::new();
        let mut res = Array1::<f64>::uninit(theta_array.len());

        for index in 0..res.len() {
            let psi = Toroidal(psi_array[index]);
            let theta = theta_array[index].rem_euclid(TAU);

            let psip = machine.qfactor().eval_other(psi, acc.xacc())?;
            let g = machine.current().eval_g(psi, acc.xacc())?;
            let b = machine.bfield().eval_b(psi, theta, acc)?;

            let kinetic = ((pzeta + psip.value()) * b / g).powi(2) / 2.0;
            let perpendicular = mu * b;
            res[index] = MaybeUninit::new(kinetic + perpendicular)
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { res.assume_init() })
    }

    /// Evaluates the Hamiltonian in a `([θ], [ψp]; Pζ, μ)` 1D array.
    ///
    /// Use this method when you *know* that `psip` is a [`Good`] coordinate.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psip` is not a [`Good`] coordinate.
    /// + Returns a [`COMError::ShapeMismatch`] if `theta_array` and `psip_array` do not have the
    ///   same length.
    /// + Returns a [`COMError::EvalError`] if `psip_array` contains out-of-bounds values.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    /// # use dexter_comspace::*;
    /// # use ndarray::Array1;
    /// # use std::f64::consts::TAU;
    /// # let path = std::path::PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, Steffen).build()?;
    /// let current = NcCurrentBuilder::new(&path, Steffen).build()?;
    /// let bfield = NcBfieldBuilder::new(&path, Bicubic).build()?;
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 10);
    /// let psip_array = Array1::linspace(0.0, machine.psip_last().value(), 10);
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let energy_array =
    ///     GCHamiltonian::eval_psip_array(machine, &theta_array, &psip_array, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_psip_array(
        machine: Machine,
        theta_array: &ArrayRef1<f64>,
        psip_array: &ArrayRef1<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array1<f64>, COMError> {
        if machine.psi_state() != FluxCoordinateState::Good {
            return Err(COMError::BadFluxCoordinate);
        }
        if !(theta_array.len() == psip_array.len()) {
            return Err(COMError::ShapeMismatch);
        }

        let mut res = Array1::<f64>::uninit(theta_array.len());

        let acc = &mut Accelerator2d::new();

        for index in 0..res.len() {
            let psip = Poloidal(psip_array[index]);
            let theta = theta_array[index].rem_euclid(TAU);

            let g = machine.current().eval_g(psip, acc.xacc())?;
            let b = machine.bfield().eval_b(psip, theta, acc)?;

            let kinetic = ((pzeta + psip.value()) * b / g).powi(2) / 2.0;
            let perpendicular = mu * b;
            res[index] = MaybeUninit::new(kinetic + perpendicular)
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { res.assume_init() })
    }

    /// Evaluates the Hamiltonian in a `([θ], [ψ]; Pζ, μ)` 2D grid.
    ///
    /// Use this method when you *know* that `psi` is a [`Good`] coordinate.
    ///
    /// The input arrays can be constructed using the [`create_poloidal_grid`] function.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psi` is not a [`Good`] coordinate.
    /// + Returns a [`COMError::ShapeMismatch`] if `theta_grid` and `psi_grid` do not have the
    ///   same shape.
    /// + Returns a [`COMError::EvalError`] if `psi_array` contains out-of-bounds values.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_comspace::*;
    /// # use ndarray::Array1;
    /// # use std::f64::consts::TAU;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 10);
    /// let psi_array = Array1::linspace(0.0, machine.psi_last().value(), 10);
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let (theta_grid, psi_grid) = create_poloidal_grid(&theta_array, &psi_array);
    /// let energy_grid =
    ///     GCHamiltonian::eval_psi_grid(machine, &theta_grid, &psi_grid, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_psi_grid(
        machine: Machine,
        theta_grid: &ArrayRef2<f64>,
        psi_grid: &ArrayRef2<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array2<f64>, COMError> {
        if machine.psi_state() != FluxCoordinateState::Good {
            return Err(COMError::BadFluxCoordinate);
        }
        if !(theta_grid.dim() == psi_grid.dim()) {
            return Err(COMError::ShapeMismatch);
        }

        let theta_array = theta_grid.row(0);
        let psi_array = psi_grid.column(0);
        let shape = (psi_array.len(), theta_array.len());
        let mut grid = Array2::<f64>::uninit(shape);

        let acc = &mut Accelerator2d::new();

        // Iterate though `psi_array` first to avoid unnecessarily recalculating `rho` and `g`.
        // Use the `Toroidal` flux for evaluation methods to avoid error propagation through double
        // interpolations.
        for m in 0..grid.nrows() {
            let psi = Toroidal(psi_array[m]);
            let psip = machine.qfactor().eval_other(psi, acc.xacc())?;
            let g = machine.current().eval_g(psi, acc.xacc())?;
            let rho = (pzeta + psip.value()) / g;
            for n in 0..grid.ncols() {
                let theta = theta_array[n].rem_euclid(TAU);
                let b = machine.bfield().eval_b(psi, theta, acc)?;
                grid[[m, n]] = MaybeUninit::new((rho * b).powi(2) / 2.0 + mu * b)
            }
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { grid.assume_init() })
    }

    /// Evaluates the Hamiltonian in a `([θ], [ψp]; Pζ, μ)` 2D grid.
    ///
    /// Use this method when you *know* that `psip` is a [`Good`] coordinate.
    ///
    /// The input arrays can be constructed using the [`create_poloidal_grid`] function.
    ///
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psip` is not a [`Good`] coordinate.
    /// + Returns a [`COMError::ShapeMismatch`] if `theta_grid` and `psip_grid` do not have the
    ///   same shape.
    /// + Returns a [`COMError::EvalError`] if `psip_array` contains out-of-bounds values.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    /// # use dexter_comspace::*;
    /// # use ndarray::Array1;
    /// # use std::f64::consts::TAU;
    /// # let path = std::path::PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, Steffen).build()?;
    /// let current = NcCurrentBuilder::new(&path, Steffen).build()?;
    /// let bfield = NcBfieldBuilder::new(&path, Bicubic).build()?;
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 10);
    /// let psip_array = Array1::linspace(0.0, machine.psip_last().value(), 10);
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let (theta_grid, psip_grid) = create_poloidal_grid(&theta_array, &psip_array);
    /// let energy_grid =
    ///     GCHamiltonian::eval_psip_grid(machine, &theta_grid, &psip_grid, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_psip_grid(
        machine: Machine,
        theta_grid: &ArrayRef2<f64>,
        psip_grid: &ArrayRef2<f64>,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array2<f64>, COMError> {
        if machine.psi_state() != FluxCoordinateState::Good {
            return Err(COMError::BadFluxCoordinate);
        }
        if !(theta_grid.dim() == psip_grid.dim()) {
            return Err(COMError::ShapeMismatch);
        }

        let theta_array = theta_grid.row(0);
        let psip_array = psip_grid.column(0);
        let shape = (psip_array.len(), theta_array.len());
        let mut grid = Array2::<f64>::uninit(shape);

        let acc = &mut Accelerator2d::new();

        for m in 0..grid.nrows() {
            let psip = Poloidal(psip_array[m]);
            let g = machine.current().eval_g(psip, acc.xacc())?;
            let rho = (pzeta + psip.value()) / g;
            for n in 0..grid.ncols() {
                let theta = theta_array[n].rem_euclid(TAU);
                let b = machine.bfield().eval_b(psip, theta, acc)?;
                grid[[m, n]] = MaybeUninit::new((rho * b).powi(2) / 2.0 + mu * b)
            }
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { grid.assume_init() })
    }
}

/// Creates a `(θ, ψ/ψp)` meshgrid from the two 1D arrays.
///
/// This function is useful for creating the meshgrid needed for the
/// [`GCHamiltonian::eval_psi_grid`] and [`GCHamiltonian::eval_psip_grid`] methods.
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
///     | E(θ0,ψ0), E(θ1,ψ0), ..., E(θM,ψ0) |
///     | E(θ0,ψ1), E(θ1,ψ1), ..., E(θM,ψ1) |
///     |     .         .     ...,     .    |
///     | E(θ0,ψN), E(θ1,ψN), ..., E(θM,ψN) |
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

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::extract::TEST_NETCDF_PATH;
    use dexter_machine::*;
    use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    use ndarray::{arr1, arr2, array};

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
        let hamiltonian_grid =
            GCHamiltonian::eval_psi_grid(machine, &theta_grid, &psi_grid, -0.03, 1e-4).unwrap();

        let expected = arr2(&[
            [0.00026296256388989685, 0.00026296256388989685],
            [0.00012897176092872728, 0.00012897176092872728],
        ])
        .to_owned();

        assert!(hamiltonian_grid.relative_eq(&expected, 1e-15, 1e-20));
    }

    #[test]
    fn toroidal_poloidal_equivalence() {
        let path = std::path::PathBuf::from(TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
        let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
        let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let pzeta = -0.0123;
        let mu = 1e-4;
        let acc = &mut Accelerator::new();
        // Avoid spline edges
        let theta_array = Array1::linspace(1e-3, TAU - 1e-3, 10);
        let psi_array = Array1::linspace(1e-3, 1.0 - 1e-3, 10) * qfactor.psi_last().value();
        let psip_array = psi_array.mapv(|psi| {
            let psi = Toroidal(psi);
            machine.qfactor().eval_other(psi, acc).unwrap().value()
        });

        let hamiltonian_psi_array =
            GCHamiltonian::eval_psi_array(machine, &theta_array, &psi_array, pzeta, mu).unwrap();
        let hamiltonian_psip_array =
            GCHamiltonian::eval_psip_array(machine, &theta_array, &psip_array, pzeta, mu).unwrap();
        assert!(hamiltonian_psi_array.relative_eq(&hamiltonian_psip_array, 1e-7, 1e-7));

        let (theta_grid, psi_grid) = create_poloidal_grid(&theta_array, &psi_array);
        let (_, psip_grid) = create_poloidal_grid(&theta_array, &psip_array);

        let hamiltonian_psi_grid =
            GCHamiltonian::eval_psi_grid(machine, &theta_grid, &psi_grid, pzeta, mu).unwrap();
        let hamiltonian_psip_grid =
            GCHamiltonian::eval_psip_grid(machine, &theta_grid, &psip_grid, pzeta, mu).unwrap();
        assert!(hamiltonian_psi_grid.relative_eq(&hamiltonian_psip_grid, 1e-7, 1e-7));
    }
}
