//! Calculations on the Constants of Motion (COMs) space.

mod energy;
mod energy_pzeta_plane;
mod tp_boundary;

use dexter_equilibrium::{Bfield, Current, FluxCommute, Qfactor};
use ndarray::{Array1, Array2};

pub use energy_pzeta_plane::EnergyPzetaPlane;
pub use tp_boundary::TrappedPassingBoundary;

use crate::COMError;

/// The constants of motion `(E, Pζ, μ)` in an unperturbed equilibrium.
#[derive(Debug)]
pub struct COMs {
    /// The Energy in Normalized Units.
    pub energy: Option<f64>,
    /// The canonical momentum `Pζ` in Normalized Units.
    pub pzeta: Option<f64>,
    /// The magnetic moment `μ` in Normalized Units.
    pub mu: Option<f64>,
}

impl COMs {
    /// Calculates the Energy on a 2D meshgrid of the `ψ` and `θ` arrays, in Normalized Units.
    ///
    /// Both `pzeta` and `mu` fields must be defined.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// # use ndarray::{Array1, Array2};
    /// # use std::f64::consts::PI;
    /// #
    /// let lcfs = LastClosedFluxSurface::Toroidal(0.1);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, lcfs);
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    ///
    /// let coms = COMs {
    ///     energy: None,
    ///     pzeta: Some(-0.025),
    ///     mu: Some(1e5),
    /// };
    ///
    /// let psi_array = Array1::linspace(0.0, lcfs.value(), 50);
    /// let theta_array = Array1::linspace(-PI, PI, 50);
    ///
    /// let grid: Array2<f64> = coms.energy_of_psi_grid(
    ///     &qfactor,
    ///     &current,
    ///     &bfield,
    ///     &psi_array,
    ///     &theta_array,
    /// )?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`COMError`] if `pzeta` or `mu` are not defined, or if any of the
    /// `psi_array` values is out of bounds.
    pub fn energy_of_psi_grid<Q, C, B>(
        &self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
        psi_array: &Array1<f64>,
        theta_array: &Array1<f64>,
    ) -> Result<Array2<f64>, COMError>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
    {
        energy::energy_of_psi_grid(self, qfactor, current, bfield, psi_array, theta_array)
    }

    /// Calculates the Energy on a 2D meshgrid of the `ψp` and `θ` arrays, in Normalized Units.
    ///
    /// Both `pzeta` and `mu` fields must be defined.
    ///
    /// Note that the [`Qfactor`] object is not needed here.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// # use std::path::PathBuf;
    /// # use ndarray::{Array1, Array2};
    /// # use std::f64::consts::PI;
    /// #
    /// let path = PathBuf::from("./netcdf.nc");
    /// let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    ///
    /// let coms = COMs {
    ///     energy: None,
    ///     pzeta: Some(-0.025),
    ///     mu: Some(1e5),
    /// };
    ///
    /// let psip_array = Array1::linspace(0.0, geometry.psip_last().unwrap(), 50);
    /// let theta_array = Array1::linspace(-PI, PI, 50);
    ///
    /// let grid: Array2<f64> = coms.energy_of_psip_grid(
    ///     &current,
    ///     &bfield,
    ///     &psip_array,
    ///     &theta_array,
    /// )?;
    /// # Ok::<_, COMError>(())
    ///
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`COMError`] if `pzeta` or `mu` are not defined, or if any of the
    /// `psi_array` values is out of bounds.
    pub fn energy_of_psip_grid<C, B>(
        &self,
        current: &C,
        bfield: &B,
        psip_array: &Array1<f64>,
        theta_array: &Array1<f64>,
    ) -> Result<Array2<f64>, COMError>
    where
        C: Current,
        B: Bfield,
    {
        energy::energy_of_psip_grid(self, current, bfield, psip_array, theta_array)
    }

    /// Constructs a [`EnergyPzetaPlane`].
    ///
    /// Only the [`COMs::mu`] field is necessary.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use dexter_simulate::*;
    /// #
    /// let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
    /// let current = LarCurrent::new();
    /// let bfield = LarBfield::new();
    ///
    /// let coms = COMs{
    ///     mu: Some(1e-5),
    ///     energy: None,
    ///     pzeta: None,
    /// };
    ///
    /// let energy_pzeta_plane = coms.build_energy_pzeta_plane(
    ///     &qfactor,
    ///     &current,
    ///     &bfield
    /// )?;
    /// # Ok::<_, COMError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`COMError`] if the [`COMs::mu`] field is `None`.
    pub fn build_energy_pzeta_plane<Q, C, B>(
        &self,
        qfactor: &Q,
        current: &C,
        bfield: &B,
    ) -> Result<EnergyPzetaPlane, COMError>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
    {
        EnergyPzetaPlane::from_coms(self, qfactor, current, bfield)
    }
}
