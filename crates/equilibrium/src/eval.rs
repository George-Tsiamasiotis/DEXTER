//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays/

use rsl_interpolation::{Accelerator, Cache};

use crate::HarmonicCache;
use crate::Result;
use crate::{Flux, Length, Radians};

/// Equilibrium geometry related quantities computation
pub trait Geometry {
    /// Calculates the radial coordinate `r(ψp)` **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let r = geom.r(0.015)?;
    /// # Ok(())
    /// # }
    /// ```
    fn r(&self, psip: Flux) -> Result<Length>;

    /// Calculates the poloidal flux `ψp(r)`, where r is **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let psip = geom.psip(0.45)?;
    /// # Ok(())
    /// # }
    /// ```
    fn psip(&self, r: Length) -> Result<Flux>;

    /// Calculates `R(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let rlab = geom.rlab(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    fn rlab(&self, psip: Flux, theta: Radians) -> Result<f64>;

    /// Calculates `Z(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let zlab = geom.zlab(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    fn zlab(&self, psip: Flux, theta: Radians) -> Result<f64>;

    /// Calculates the Jacobian `J(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let geom = Geometry::from_dataset(&path, "steffen", "bicubic")?;
    ///
    /// let j = geom.jacobian(0.015, 2.0*PI)?;
    /// # Ok(())
    /// # }
    /// ```
    fn jacobian(&self, psip: Flux, theta: Radians) -> Result<f64>;
}

/// q-factor related quantities computation.
pub trait Qfactor {
    /// Calculates the q-factor `q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let q = qfactor.q(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn q(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the toroidal flux `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let psi = qfactor.psi(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn psi(&self, psip: Flux, acc: &mut Accelerator) -> Result<Flux>;

    /// Calculates the derivative `dψ/dψp`.
    ///
    /// This value should always equal `q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let qfactor = Qfactor::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dpsi_dpsip = qfactor.dpsi_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dpsi_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Plasma current related quantities computation.
pub trait Current {
    /// Calculates `g(ψp)`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let g = currents.g(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn g(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `I(ψp)`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let i = currents.i(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn i(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕g(ψp)/𝜕ψp`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let dg = currents.dg_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dg_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕I(ψp)/𝜕ψp`
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let currents = Currents::from_dataset(&path, "cubic")?;
    ///
    /// let mut acc = Accelerator::new();
    /// let di = currents.di_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn di_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Magnetic field related quantities computation.
pub trait Bfield {
    /// Calculates `B(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let b =  bfield.b(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn b(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dpsip = bfield.db_dpsip(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn db_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let db_dtheta = bfield.db_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn db_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕²B(ψp, θ) /𝜕𝜓p²`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn d2b_dpsip2(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕²B(ψp, θ) /𝜕θ²`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dtheta2(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn d2b_dtheta2(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕²B(ψp, θ) /𝜕ψp𝜕θ`.
    ///
    /// # Note:
    ///
    /// The data are calculated by evaluating the bfield spline's derivative, rather than
    /// extracting the data arrays from the netCDF file.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let bfield = Bfield::from_dataset(&path, "bicubic")?;
    ///
    /// let mut psi_acc = Accelerator::new();
    /// let mut theta_acc = Accelerator::new();
    /// let mut cache = Cache::new();
    ///
    /// let d2b_dpsip2 = bfield.d2b_dpsip_dtheta(0.015, 2.0*PI, &mut psi_acc, &mut theta_acc, &mut cache)?;
    /// # Ok(())
    /// # }
    /// ```
    fn d2b_dpsip_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;
}

/// Single Harmonic related quantities computation
pub trait Harmonic {
    /// Calculates the harmonic `α(ψp) * cos(mθ-nζ+φ(ψp))`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let h = harmonic.h(0.015, 2.0*PI, 0.0, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn h(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dpsip = harmonic.dh_dpsip(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dh_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dtheta = harmonic.dh_dtheta(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dh_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕ζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dzeta = harmonic.dh_dzeta(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dh_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕t`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dt = harmonic.dh_dt(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    #[allow(unused_variables)]
    fn dh_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        // Time-independent perturbations at the moment.
        Ok(0.0)
    }

    /// Calculates the harmonic's *amplitude* `α(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let a = harmonic.a(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn a(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *amplitude* derivative `dα(ψp)/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let da_dpsip = harmonic.da_dpsip(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn da_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *phase* `φ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let phase = harmonic.phase(0.015, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn phase(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64>;
}

/// Perturbation related quantities computation
pub trait Perturbation {
    /// Calculates the Perturbation `Σ{ α(n,m)(ψp) * cos(mθ-nζ+φ0) }`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let p = per.p(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn p(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ψp`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dpsip = per.dp_dpsip(0.015, 2.0*PI, PI,&mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dp_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dtheta = per.dp_dtheta(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dp_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dzeta = per.dp_dzeta(0.015, 2.0*PI, PI, &mut hcaches, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dp_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `t`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::*;
    /// # use std::f64::consts::PI;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonics = vec![
    ///     Harmonic::from_dataset(&path, "akima", 3, 1)?,
    ///     Harmonic::from_dataset(&path, "akima", 3, 2)?,
    /// ];
    /// let per = Perturbation::from_harmonics(&harmonics);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcache = vec![HarmonicCache::new(); harmonics.len()];
    /// let dp_dt = per.dp_dt(0.015, 2.0*PI, PI, &mut hcache, &mut acc)?;
    /// # Ok(())
    /// # }
    /// ```
    fn dp_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        caches: &mut [Box<dyn HarmonicCache>],
        acc: &mut Accelerator,
    ) -> Result<f64>;
}
