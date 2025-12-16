//! Definitions of evaluation methods of equilibrium objects.
//!
//! For analytical equilibria, this is achieved by evaluation of analytical formulas, while for
//! numerical equilibria by interpolation over the reconstructed data arrays.

use rsl_interpolation::{Accelerator, Cache};

use crate::Result;
use crate::cache::*;

// TODO: (maybe) add doctests

/// Equilibrium geometry related quantities computation
pub trait Geometry {
    /// Calculates the radial coordinate `r(ψp)` **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let r = geometry.r(0.015)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn r(&self, psip: f64) -> Result<f64>;

    /// Calculates the poloidal flux `ψp(r)`, where r is **in \[m\]**.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let psip = geometry.psip(0.02)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psip(&self, r: f64) -> Result<f64>;

    /// Calculates the toroidal flux `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// let psi = geometry.psi(0.015)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psi(&self, psip: f64) -> Result<f64>;

    /// Calculates `R(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let rlab = geometry.rlab(0.015, 3.1415)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn rlab(&self, psip: f64, theta: f64) -> Result<f64>;

    /// Calculates `Z(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let zlab = geometry.zlab(0.015, 3.1415)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn zlab(&self, psip: f64, theta: f64) -> Result<f64>;

    /// Calculates the Jacobian `J(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let geometry = NcGeometryBuilder::new(&path, "steffen", "bicubic").build()?;
    /// #
    /// let j = geometry.jacobian(0.015, 3.1415)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn jacobian(&self, psip: f64, theta: f64) -> Result<f64>;
}

/// q-factor related quantities computation.
pub trait Qfactor: Clone {
    /// Calculates the q-factor `q(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let q = qfactor.q(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn q(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the toroidal flux `ψ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let psi = qfactor.psi(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn psi(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the derivative `dψ/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let qfactor = NcQfactorBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let dpsi_dpsip = qfactor.dpsi_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;
}

/// Plasma current related quantities computation.
pub trait Current {
    /// Calculates `g(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let g = current.g(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn g(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `I(ψp)`
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let current = NcCurrentBuilder::new(&path, "steffen").build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let i = current.i(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn i(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates `𝜕g(ψp)/𝜕ψp`
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
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

    /// Calculates `𝜕I(ψp)/𝜕ψp`
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
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

/// Magnetic field related quantities computation.
pub trait Bfield {
    /// Calculates `B(ψp, θ)`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bilinear").build()?;
    /// #
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let b = bfield.b(0.015, 3.1415, &mut xacc, &mut yacc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn b(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bilinear").build()?;
    /// #
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dpsip = bfield.db_dpsip(0.015, 3.1415, &mut xacc, &mut yacc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_dpsip(
        &self,
        psip: f64,
        theta: f64,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
        cache: &mut Cache<f64>,
    ) -> Result<f64>;

    /// Calculates `𝜕B(ψp, θ) /𝜕𝜃`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, "bilinear").build()?;
    /// #
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    /// let mut cache = Cache::new();
    /// let db_dtheta = bfield.db_dtheta(0.015, 3.1415, &mut xacc, &mut yacc, &mut cache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn db_dtheta(
        &self,
        psip: f64,
        theta: f64,
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
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let h = harmonic.h(0.015, 3.1415, 6.2831, &mut acc, &mut hcache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn h(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕ψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dpsip = harmonic.dh_dpsip(0.015, 3.1415, 6.2831, &mut acc, &mut hcache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dh_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64>;

    /// Calculates the harmonic derivative `𝜕h/𝜕θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dtheta = harmonic.dh_dtheta(0.015, 3.1415, 6.2831, &mut acc, &mut hcache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dh_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕ζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dzeta = harmonic.dh_dzeta(0.015, 3.1415, 6.2831, &mut acc, &mut hcache)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dh_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64>;

    /// Calculates the perturbation derivative `𝜕h/𝜕t`.
    ///
    /// Time-independent perturbations at the moment, so it always returns `0.0`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let mut hcache = HarmonicCache::new();
    /// let dh_dt = harmonic.dh_dt(0.015, 3.1415, 6.2831, &mut acc, &mut hcache)?;
    /// assert_eq!(dh_dt, 0.0);
    /// # Ok::<_, EqError>(())
    /// ```
    #[allow(unused_variables)]
    fn dh_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64> {
        Ok(0.0)
    }

    /// Calculates the harmonic's *amplitude* `α(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let a = harmonic.a(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn a(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *amplitude* derivative `dα(ψp)/dψp`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let da_dpsip = harmonic.da_dpsip(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn da_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the harmonic's *phase* `φ(ψp)`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let phase = harmonic.phase(0.015, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn phase(&self, psip: f64, acc: &mut Accelerator) -> Result<f64>;

    /// Calculates the term inside the cosine, modulo 2π.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?;
    /// #
    /// let mut acc = Accelerator::new();
    /// let module = harmonic.mod_arg(0.015, 3.1415, 6.2831, &mut acc)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn mod_arg(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64>;
}

/// Perturbation related quantities computation
///
/// Since a Perturbation consists of multiple harmonics, the perturbation's value at a specific
/// point is the total sum of all the contained harmonics evaluated at that point.
pub trait Perturbation {
    /// Calculates the perturbation's value.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = NcPerturbation::from_harmonics(&vec![
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    /// ]);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = [HarmonicCache::new(), HarmonicCache::new()];
    /// let p = perturbation.p(0.015, 3.1415, 6.2831, &mut acc, &mut hcaches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn p(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ψp`,
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = NcPerturbation::from_harmonics(&vec![
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    /// ]);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = [HarmonicCache::new(), HarmonicCache::new()];
    /// let dp_dpsip = perturbation.dp_dpsip(0.015, 3.1415, 6.2831, &mut acc, &mut hcaches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dp_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `θ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = NcPerturbation::from_harmonics(&vec![
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    /// ]);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = [HarmonicCache::new(), HarmonicCache::new()];
    /// let dp_dtheta = perturbation.dp_dtheta(0.015, 3.1415, 6.2831, &mut acc, &mut hcaches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dp_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `ζ`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = NcPerturbation::from_harmonics(&vec![
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    /// ]);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = [HarmonicCache::new(), HarmonicCache::new()];
    /// let dp_dzeta = perturbation.dp_dzeta(0.015, 3.1415, 6.2831, &mut acc, &mut hcaches)?;
    /// # Ok::<_, EqError>(())
    /// ```
    fn dp_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> Result<f64>;

    /// Calculates the Perturbation's derivative with respect to `t`.
    ///
    /// Time-independent perturbations at the moment, so it always returns `0.0`.
    ///
    /// # Example
    ///
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::{Accelerator, Cache};
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let perturbation = NcPerturbation::from_harmonics(&vec![
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    ///    NcHarmonicBuilder::new(&path, "steffen", 1, 2).build()?,
    /// ]);
    ///
    /// let mut acc = Accelerator::new();
    /// let mut hcaches = [HarmonicCache::new(), HarmonicCache::new()];
    /// let dp_dt = perturbation.dp_dt(0.015, 3.1415, 6.2831, &mut acc, &mut hcaches)?;
    /// assert_eq!(dp_dt, 0.0);
    /// # Ok::<_, EqError>(())
    /// ```
    fn dp_dt(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        caches: &mut [HarmonicCache],
    ) -> Result<f64>;

    /// Returns the number of harmonics.
    fn len(&self) -> usize;

    /// Returns true if the perturbation has no harmonics (== no perturbation).
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
