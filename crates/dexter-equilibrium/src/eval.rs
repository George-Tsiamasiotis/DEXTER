//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays.

use ndarray::Array1;
use rsl_interpolation::{Accelerator, Cache};

use crate::Result;

/// Equilibrium geometry related quantities computation.
pub trait Geometry {
    /// Calculates the radial coordinate `r(ψ)` **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let r_of_psi = geometry.r_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn r_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the radial coordinate `r(ψp)` **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let r_of_psip = geometry.r_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn r_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `ψ(r)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi_of_r = geometry.psi_of_r(0.02, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psi_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `ψp(r)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psip_of_r = geometry.psip_of_r(0.02, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `R(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let rlab_of_psi = geometry.rlab_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn rlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `R(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let rlab_of_psip = geometry.rlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn rlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `Z(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let zlab_of_psi = geometry.zlab_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn zlab_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `Z(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let zlab_of_psip = geometry.zlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn zlab_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates the Jacobian `J(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let jacobian_of_psi = geometry.jacobian_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn jacobian_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates the Jacobian `J(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut psip_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let jacobian_of_psip = geometry.zlab_of_psip(0.01, 3.14, &mut psip_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn jacobian_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psip_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Returns the last `Rlab` values that correspond to the device's wall.
    fn rlab_wall(&self) -> Array1<f64>;

    /// Returns the last `Zlab` values that correspond to the device's wall.
    fn zlab_wall(&self) -> Array1<f64>;
}

/// Conversion between the two flux coordinates `ψ` and `ψp`
pub trait FluxCommute {
    /// Calculates `ψp(ψ)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psip_of_psi = qfactor.psip_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi_of_psip = geometry.psi_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;
}

/// Magnetic field related quantities computation
pub trait Bfield {
    /// Calculates `B(ψ, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let b_of_psi = bfield.b_of_psi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn b_of_psi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `B(ψp, θ)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let b_of_psip = bfield.b_of_psip(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn b_of_psip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `dB(ψ, θ)/dψ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dpsi = bfield.db_dpsi(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_dpsi(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `dB(ψp, θ)/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dpsip = bfield.db_dpsip(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `dB(ψ, θ)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_of_psi_dtheta = bfield.db_of_psi_dtheta(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `dB(ψp, θ)/dθ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bicubic").build()?;
    /// #
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_of_psip_dtheta = bfield.db_of_psip_dtheta(0.01, 3.14, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        psi_acc: &mut Accelerator,
        theta_acc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;
}

/// q-factor related quantities computation.
pub trait Qfactor {
    /// Calculates `q(ψ)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let q_of_psi = qfactor.q_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `q(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let q_of_psip = qfactor.q_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the derivative `dψp(ψ)/dψ`.
    ///
    /// It's a good check that the values coincide with `qfactor.iota_of_psi(psi)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dpsip_dpsi = qfactor.dpsip_dpsi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the derivative `dψ(ψp)/dψp`.
    ///
    /// It's a good check that the values coincide with `qfactor.q_of_psip(psip)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dpsi_dpsip = qfactor.dpsi_dpsip(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `i(ψ) = 1/q(ψ)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let iota_of_psi = qfactor.iota_of_psi(0.01, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn iota_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.q_of_psi(psi, acc)?.recip())
    }

    /// Calculates `i(ψp) = 1/q(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let iota_of_psip = qfactor.iota_of_psip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn iota_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.q_of_psip(psip, acc)?.recip())
    }
}

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
