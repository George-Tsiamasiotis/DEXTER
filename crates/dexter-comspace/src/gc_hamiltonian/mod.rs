//! Defines the unit struct [`GCHamiltonian`] and several basic methods.
//!
//! `impl` blocks are separated and grouped into several submodules for clarity.

#![expect(clippy::min_ident_chars, reason = "Hamiltonian terms")]

mod array;
mod curve;
mod point;
mod segment;

use std::f64::consts::TAU;

pub use curve::StationaryCurve;
pub use point::{CriticalPoint, CriticalPointKind};
pub use segment::ParametrizedSegment;

// ===============================================================================================

use dexter_machine::{Accelerator2d, EvalError, Machine, MagneticFlux};

/// Guiding center Hamiltonian.
///
/// Provides evaluation methods of the Hamiltonian and its derived quantities.
///
/// This unit struct has no internal state; the provided methods require all parameters to be passed
/// as arguments and calculate all needed quantities in every call. Consider implementing optimized
/// versions of these calculations in performance-critical routines.
#[derive(Debug, Clone, Copy)]
pub struct GCHamiltonian;

impl GCHamiltonian {
    /// Evaluates the Hamiltonian in a `(θ, ψ/ψp; Pζ, μ)` point.
    ///
    /// # Errors
    ///
    /// + This method returns an [`EvalError`] if `flux` is out of bounds.
    /// + This method returns an [`EvalError::UndefinedEvaluation`] if a flux conversion is need but
    ///   not defined.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use dexter_comspace::*;
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.05);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    ///
    /// let theta = 3.1415;
    /// let psi = MagneticFlux::Toroidal(0.015);
    /// let pzeta = -0.1*qfactor.psip_last().value();
    /// let mu = 1e-4;
    ///
    /// let energy: f64 = GCHamiltonian::eval(machine, theta, psi, pzeta, mu)?;
    /// # Ok::<_, EvalError>(())
    ///
    /// ```
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn eval(
        machine: Machine,
        theta: f64,
        flux: MagneticFlux,
        pzeta: f64,
        mu: f64,
    ) -> Result<f64, EvalError> {
        use dexter_machine::MagneticFluxKind::*;

        let acc = &mut Accelerator2d::new();
        let theta = theta.rem_euclid(TAU);

        let (b, g, psip): (f64, f64, MagneticFlux);
        match (machine.good_flux(), flux.kind()) {
            // Use the passed `flux` in the `b` and `g` calculations to avoid error propagation
            (Toroidal, Toroidal) => {
                b = machine.bfield().eval_b(flux, theta, acc)?;
                g = machine.current().eval_g(flux, acc.xacc())?;
                psip = machine.qfactor().eval_other(flux, acc.xacc())?;
            }
            (Poloidal, Poloidal) => {
                b = machine.bfield().eval_b(flux, theta, acc)?;
                g = machine.current().eval_g(flux, acc.xacc())?;
                psip = flux;
            }
            // Convert to the other flux and use `machine.good_flux()` in the `b` and `g` calculations
            (Toroidal, Poloidal) => {
                let psi = machine.qfactor().eval_other(flux, acc.xacc())?;
                b = machine.bfield().eval_b(psi, theta, acc)?;
                g = machine.current().eval_g(psi, acc.xacc())?;
                psip = flux;
            }
            (Poloidal, Toroidal) => {
                psip = machine.qfactor().eval_other(flux, acc.xacc())?;
                b = machine.bfield().eval_b(psip, theta, acc)?;
                g = machine.current().eval_g(psip, acc.xacc())?;
            }
        };

        let kinetic = ((pzeta + psip.value()) * b / g).powi(2) / 2.0;
        let perpendicular = mu * b;

        Ok(kinetic + perpendicular)
    }

    /// Builds the [`StationaryCurve`] of `machine`.
    ///
    /// The `StationaryCurve` is constructed w.r.t [`Machine::good_flux`].
    #[must_use]
    pub fn build_stationary_curve(machine: Machine) -> StationaryCurve {
        StationaryCurve::build(machine)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::MagneticFlux::*;
    use dexter_machine::extract::{POLOIDAL_TEST_NETCDF_PATH, TOROIDAL_TEST_NETCDF_PATH};
    use dexter_machine::*;
    use dexter_machine::{Interpolation1dType::Steffen, Interpolation2dType::Bicubic};

    #[test]
    fn gc_hamiltonian_evaluate_machine_good_both() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let _ = GCHamiltonian::eval(machine, 0.0, Toroidal(0.01), -0.01, 1e-4).unwrap();
        let _ = GCHamiltonian::eval(machine, 0.0, Poloidal(0.01), -0.01, 1e-4).unwrap();
    }

    #[test]
    fn gc_hamiltonian_evaluate_machine_good_psi() {
        let path = std::path::PathBuf::from(TOROIDAL_TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
        let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
        let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let _ = GCHamiltonian::eval(machine, 0.0, Toroidal(0.01), -0.01, 1e-4).unwrap();
        // `psip` is `Bad` so we cannot calculate `psi`
        let _ = GCHamiltonian::eval(machine, 0.0, Poloidal(0.01), -0.01, 1e-4).unwrap_err();
    }

    #[test]
    fn gc_hamiltonian_evaluate_machine_good_psip() {
        let path = std::path::PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let qfactor = NcQfactorBuilder::new(&path, Steffen).build().unwrap();
        let current = NcCurrentBuilder::new(&path, Steffen).build().unwrap();
        let bfield = NcBfieldBuilder::new(&path, Bicubic).build().unwrap();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let _ = GCHamiltonian::eval(machine, 0.0, Poloidal(0.01), -0.01, 1e-4).unwrap();
        // `psi` is `Bad` so we cannot calculate `psip`
        let _ = GCHamiltonian::eval(machine, 0.0, Toroidal(0.01), -0.01, 1e-4).unwrap_err();
    }
}
