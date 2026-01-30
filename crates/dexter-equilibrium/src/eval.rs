//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays.

use rsl_interpolation::Accelerator;

use crate::Result;

/// Plasma current related quantities computation.
pub trait Current {
    /// Calculates `g(ψ)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let g_of_psi = current.g_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn g_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `g(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let g_of_psip = current.g_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn g_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `I(ψ)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let i_of_psi = current.i_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn i_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `I(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let i_of_psip = current.i_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn i_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕g(ψ)/𝜕ψ`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dg_dpsi = current.dg_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dg_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕g(ψp)/𝜕ψp`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dg_dpsip = current.dg_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dg_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕I(ψ)/𝜕ψ`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let di_dpsi = current.di_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn di_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕I(ψp)/𝜕ψp`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let di_dpsip = current.di_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn di_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;
}
