//! Evaluation of the Hamiltonian in 1D arrays and 2D meshgrids.

#![expect(clippy::min_ident_chars, reason = "Hamiltonian terms")]

use ndarray::{Array1, Array2, ArrayRef1};
use rsl_interpolation::Accelerator2d;
use std::{f64::consts::TAU, mem::MaybeUninit};

use dexter_machine::{FluxCoordinateState, Machine, MagneticFluxArray, MagneticFluxKind};

use crate::{COMError, GCHamiltonian};

impl GCHamiltonian {
    /// Evaluates the Hamiltonian in a `([θ], [ψ/ψp]; Pζ, μ)` 1D array.
    ///
    /// The produced array has the form:
    ///
    /// ```text
    /// | H(θ0,ψ0), H(θ1,ψ1), ..., H(θN,ψN) |
    /// ```
    ///
    /// This method uses the [`flux_array.kind()`][MagneticFluxKind] flux coordinate to perform all
    /// evaluations. Since only `ψp` appears explicitly in the Hamiltonian, if
    /// `flux_array.kind() == Toroidal` then an extra evaluation is needed to calculate `ψp`.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psi` is not a [`Good`] coordinate of `machine`.
    /// + Returns a [`COMError::ShapeMismatch`] if `theta_array` and `flux_array` do not have the
    ///   same length.
    /// + Returns a [`COMError::EvalError`] if `flux_array` contains any out-of-bounds values or
    ///   if any evaluation fails.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_comspace::*;
    /// # use ndarray::*;
    /// # use std::f64::consts::TAU;
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.9, MagneticFlux::Toroidal(0.05));
    /// # let current = LarCurrent::new();
    /// # let bfield = LarBfield::new();
    /// # let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    /// #
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 10);
    /// let psi_array = MagneticFluxArray::Toroidal(
    ///     Array1::linspace(1e-3, 1.0 - 1e-3, 10) * qfactor.psi_last().value(),
    /// );
    ///
    /// let energy_array: Array1<f64> =
    ///     GCHamiltonian::eval_array(machine, &theta_array, &psi_array, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_array(
        machine: Machine,
        theta_array: &ArrayRef1<f64>,
        flux_array: &MagneticFluxArray,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array1<f64>, COMError> {
        if !machine.is_good_flux(flux_array.kind()) {
            return Err(COMError::BadFluxCoordinate);
        }
        if !(theta_array.len() == flux_array.len()) {
            return Err(COMError::ShapeMismatch);
        }

        let acc = &mut Accelerator2d::new();
        let mut array = Array1::<f64>::uninit(theta_array.len());

        // Use the passed `flux_array` to evaluate all machine quantities, and only calculate `psip`
        // if it is necessary, that is, `flux_array.kind() == Toroidal`
        for index in 0..array.len() {
            let theta = theta_array[index].rem_euclid(TAU);
            let Some(flux) = flux_array.get(index) else {
                unreachable!("index is in-bounds")
            };
            let psip = match flux.kind() {
                MagneticFluxKind::Toroidal => flux.to_other(machine.qfactor())?,
                MagneticFluxKind::Poloidal => flux,
            };

            let g = machine.current().eval_g(flux, acc.xacc())?;
            let b = machine.bfield().eval_b(flux, theta, acc)?;

            let kinetic = ((pzeta + psip.value()) * b / g).powi(2) / 2.0;
            let perpendicular = mu * b;
            array[index] = MaybeUninit::new(kinetic + perpendicular)
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { array.assume_init() })
    }

    /// Evaluates the Hamiltonian in a `([θ], [ψ/ψp]; Pζ, μ)` 2D grid.
    ///
    /// The `theta` and `flux` arrays do not need to have the same length. The produced grid has the form:
    ///
    /// ```text
    /// | H(θ0,ψ0), H(θ1,ψ0), ..., H(θM,ψ0) |
    /// | H(θ0,ψ1), H(θ1,ψ1), ..., H(θM,ψ1) |
    /// |     .         .     ...,     .    |
    /// | H(θ0,ψN), H(θ1,ψN), ..., H(θM,ψN) |
    /// ```
    ///
    /// This method uses the [`flux_array.kind()`][MagneticFluxKind] flux coordinate to perform all
    /// evaluations. Since only `ψp` appears explicitly in the Hamiltonian, if
    /// `flux_array.kind() == Toroidal` then an extra evaluation is needed to calculate `ψp`.
    ///
    /// The only allocation that takes place is the allocation of the returned grid. No extra
    /// allocations for the `(θ, ψ/ψp)` meshgrid are needed.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `psi` is not a [`Good`] coordinate of `machine`.
    /// + Returns a [`COMError::EvalError`] if `flux_array` contains any out-of-bounds values or
    ///   if any evaluation fails.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_comspace::*;
    /// # use ndarray::*;
    /// # use std::f64::consts::TAU;
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.9, MagneticFlux::Toroidal(0.05));
    /// # let current = LarCurrent::new();
    /// # let bfield = LarBfield::new();
    /// # let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    /// #
    /// let pzeta = -0.1 * machine.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let theta_array = Array1::linspace(0.0, TAU, 20);
    /// let psi_array = MagneticFluxArray::Toroidal(
    ///     Array1::linspace(1e-3, 1.0 - 1e-3, 10) * qfactor.psi_last().value(),
    /// );
    ///
    /// let energy_grid: Array2<f64> =
    ///     GCHamiltonian::eval_grid(machine, &theta_array, &psi_array, pzeta, mu)?;
    /// # Ok::<_, COMError>(())
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval_grid(
        machine: Machine,
        theta_array: &ArrayRef1<f64>,
        flux_array: &MagneticFluxArray,
        pzeta: f64,
        mu: f64,
    ) -> Result<Array2<f64>, COMError> {
        if machine.psi_state() != FluxCoordinateState::Good {
            return Err(COMError::BadFluxCoordinate);
        }

        let shape = (flux_array.len(), theta_array.len());
        let mut grid = Array2::<f64>::uninit(shape);

        let acc = &mut Accelerator2d::new();

        // Use the passed `flux_array` to evaluate all machine quantities, and only calculate `psip`
        // if it is necessary, that is, `flux_array.kind() == Toroidal`
        //
        // Iterate through the rows first: `rho` only depends on `flux`, therefore we only need to
        // calculate it once per row.
        for m in 0..grid.nrows() {
            let Some(flux) = flux_array.get(m) else {
                unreachable!("index is in-bounds")
            };
            let psip = match flux.kind() {
                MagneticFluxKind::Toroidal => flux.to_other(machine.qfactor())?,
                MagneticFluxKind::Poloidal => flux,
            };

            let g = machine.current().eval_g(flux, acc.xacc())?;
            let rho = (pzeta + psip.value()) / g;
            for n in 0..grid.ncols() {
                let theta = theta_array[n].rem_euclid(TAU);
                let b = machine.bfield().eval_b(flux, theta, acc)?;
                grid[[m, n]] = MaybeUninit::new((rho * b).powi(2) / 2.0 + mu * b)
            }
        }

        // SAFETY: The loop passes from all elements and initializes them
        Ok(unsafe { grid.assume_init() })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::extract::TEST_NETCDF_PATH;
    use dexter_machine::*;
    use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};
    use ndarray::{arr1, arr2, array};

    #[test]
    fn gcmotion_check() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let theta_array = arr1(&vec![-1.0, 1.0]);
        let psi_array = MagneticFluxArray::Toroidal(array![0.01, 0.02]);

        let hamiltonian_grid =
            GCHamiltonian::eval_grid(machine, &theta_array, &psi_array, -0.03, 1e-4).unwrap();

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
        // Avoid spline edges
        let theta_array = Array1::linspace(1e-3, TAU - 1e-3, 10);
        let psi_array = MagneticFluxArray::Toroidal(
            Array1::linspace(1e-3, 1.0 - 1e-3, 10) * qfactor.psi_last().value(),
        );
        let psip_array = psi_array.to_other(&qfactor).unwrap();

        let hamiltonian_psi_array =
            GCHamiltonian::eval_array(machine, &theta_array, &psi_array, pzeta, mu).unwrap();
        let hamiltonian_psip_array =
            GCHamiltonian::eval_array(machine, &theta_array, &psip_array, pzeta, mu).unwrap();
        assert!(hamiltonian_psi_array.relative_eq(&hamiltonian_psip_array, 1e-7, 1e-7));

        let psi_array = MagneticFluxArray::Toroidal(
            Array1::linspace(1e-3, 1.0 - 1e-3, 20) * qfactor.psi_last().value(),
        );
        let psip_array = psi_array.to_other(&qfactor).unwrap();

        let hamiltonian_psi_grid =
            GCHamiltonian::eval_grid(machine, &theta_array, &psi_array, pzeta, mu).unwrap();
        let hamiltonian_psip_grid =
            GCHamiltonian::eval_grid(machine, &theta_array, &psip_array, pzeta, mu).unwrap();
        assert!(hamiltonian_psi_grid.relative_eq(&hamiltonian_psip_grid, 1e-7, 1e-7));
    }
}
