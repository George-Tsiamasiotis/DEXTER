use std::f64::consts::TAU;
use std::path::PathBuf;

use rsl_interpolation::{Accelerator, DynSpline, make_spline};
use utils::array1D_getter_impl;

use crate::Result;
use crate::{Flux, Length, Radians};

use ndarray::Array1;

/// Single perturbation harmonic reconstructed from a netCDF file.
///
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ(ψp))`, where `α(ψp)` and `φ(ψp)` are calculated by
/// interpolation over numerical data.
pub struct Harmonic {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,
    /// Spline over the perturbation amplitude `α` data, as a function of ψp.
    a_spline: DynSpline<f64>,
    /// Spline over the perturbation amplitude `φ` data, as a function of ψp.
    phase_spline: DynSpline<f64>,
    /// The mean value of the phase data array.
    ///
    /// This value is used when the [`phase-interpolation`] feature is enabled.
    ///
    /// [`phase-interpolation`]: index.html#features
    phase_average: Radians,
    /// The `θ` frequency number.
    m: i64,
    /// The `ζ` frequency number.
    n: i64,

    // Used in the actual calculations
    _m: f64,
    _n: f64,
}

/// Holds the Harmonic's values evalutated at a specific point.
///
/// Since all the harmonic's methods are called consecutively over the same coordinates, most terms
/// do not need to be calculated every time.
///
/// Similar to the Accelerators, they are stored inside State, and do not affect the behavior of the
/// equilibrium objects themselves.
///
/// The cache should be cloned in each new state calculated from the Solver.
#[derive(Clone, Default)]
pub struct HarmonicCache {
    hits: usize,
    misses: usize,
    psip: Flux,
    theta: Radians,
    zeta: Radians,
    pub(crate) alpha: f64,
    pub(crate) phase: Radians,
    pub(crate) dalpha: f64,
    pub(crate) sin: f64,
    pub(crate) cos: f64,
}

impl HarmonicCache {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn hits(&self) -> usize {
        self.hits
    }

    pub fn misses(&self) -> usize {
        self.misses
    }

    /// Checks if the cache's fields are valid.
    ///
    /// Comparing floats is OK here since they are simply copied between every call, and we want
    /// the check to fail with the slightest difference.
    pub(crate) fn is_updated(&mut self, psip: Flux, theta: Radians, zeta: Radians) -> bool {
        if (self.psip == psip) & (self.theta == theta) & (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    /// Updates the cache's fields.
    pub(crate) fn update(
        &mut self,
        h: &Harmonic,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        acc: &mut Accelerator,
    ) -> Result<()> {
        self.psip = psip;
        self.theta = theta;
        self.zeta = zeta;
        self.alpha = h.a_spline.eval(psip, acc)?;
        self.phase = calculate_phase(h, psip, acc)?;
        self.dalpha = h.a_spline.eval_deriv(psip, acc)?;
        let mod_arg = (h._m * self.theta - h._n * self.zeta + self.phase).rem_euclid(TAU);
        (self.sin, self.cos) = mod_arg.sin_cos();
        Ok(())
    }
}

/// Returns the phase by interpolating over the extracted phase data.
#[cfg(feature = "phase-interpolation")]
#[inline(always)]
fn calculate_phase(h: &Harmonic, psip: f64, acc: &mut Accelerator) -> Result<Radians> {
    Ok(h.phase_spline.eval(psip, acc)?)
}

/// Simply returns the harmonic's average phase.
#[cfg(not(feature = "phase-interpolation"))]
#[allow(unused_variables)]
#[inline(always)]
fn calculate_phase(h: &Harmonic, psip: f64, acc: &mut Accelerator) -> Result<Radians> {
    Ok(h.phase_average)
}

// Creation
impl Harmonic {
    /// Constructs a [`Harmonic`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// The spline is only over the amplitude `α` and phase `φ`, of the perturbation. The rest of the
    /// exrpession is analytic.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = Harmonic::from_dataset(&path, "akima", 3, 2)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_dataset(path: &PathBuf, typ: &str, m: i64, n: i64) -> Result<Self> {
        use crate::extract::*;
        use config::netcdf_fields::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_owned();
        let (a_data, phase_data) = extract_harmonic_arrays(&f, m, n)?;

        // `extract_array()` has already checked if the arrays are empty
        let a_spline = make_spline(
            typ,
            psip_data.as_slice().expect("array is non-empty"),
            a_data.as_slice().expect("array is non-empty"),
        )?;
        // We still want the phase spline for plotting, even with the 'phase-average' feature
        // enabled.
        let phase_spline = make_spline(
            typ,
            psip_data.as_slice().expect("array is non-empty"),
            phase_data.as_slice().expect("array is non-empty"),
        )?;
        let phase_average = phase_data.mean().expect("array is non-empty");

        Ok(Self {
            path: path.to_owned(),
            typ: typ.into(),
            a_spline,
            phase_spline,
            phase_average,
            m,
            n,
            _m: m as f64,
            _n: n as f64,
        })
    }
}

// Interpolation
impl Harmonic {
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
    pub fn h(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * cache.cos)
    }

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
    pub fn dh_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.dalpha * cache.cos)
    }

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
    pub fn dh_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * (-self._m) * cache.sin)
    }

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
    pub fn dh_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * self._n * cache.sin)
    }

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
    pub fn dh_dt(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut HarmonicCache,
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
    pub fn a(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.a_spline.eval(psip, acc)?)
    }

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
    pub fn da_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.a_spline.eval_deriv(psip, acc)?)
    }

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
    pub fn phase(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self.phase_spline.eval(psip, acc)?)
    }
}

/// Getters
impl Harmonic {
    /// Returns the netCDF file's path.
    pub fn path(&self) -> PathBuf {
        self.path.clone()
    }

    /// Returns the interpolation type.
    pub fn typ(&self) -> String {
        self.typ.clone()
    }

    /// Returns the number of data points.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.a_spline.xa.len()
    }

    /// Returns the poloidal mode number `m`.
    pub fn m(&self) -> i64 {
        self.m
    }

    /// Returns the poloidal mode number `m`.
    pub fn n(&self) -> i64 {
        self.n
    }

    /// Returns the phase average.
    pub fn phase_average(&self) -> f64 {
        self.phase_average
    }

    array1D_getter_impl!(psip_data, a_spline.xa, Flux);
    array1D_getter_impl!(a_data, a_spline.ya, Length);
    array1D_getter_impl!(phase_data, phase_spline.ya, Radians);
}

impl Clone for Harmonic {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            a_spline: make_spline(&self.typ, &self.a_spline.xa, &self.a_spline.ya)
                .expect("Could not clone spline."),
            phase_spline: make_spline(&self.typ, &self.phase_spline.xa, &self.phase_spline.ya)
                .expect("Could not clone spline."),
            phase_average: self.phase_average,
            m: self.m,
            n: self.n,
            _m: self._m,
            _n: self._n,
        }
    }
}

impl std::fmt::Debug for Harmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Harmonic")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("m", &self.m())
            .field("n", &self.n())
            .field("phase_average", &format!("{:.7}", self.phase_average()))
            .finish()
    }
}
impl std::fmt::Debug for HarmonicCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HarmonicCache")
            .field("hits  ", &self.hits)
            .field("misses", &self.misses)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use config::STUB_NETCDF_PATH;

    fn create_harmonic() -> Harmonic {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        Harmonic::from_dataset(&path, "steffen", 0, 1).unwrap()
    }

    #[test]
    fn test_creation() {
        let h = create_harmonic();
        let _ = h.clone();
        let _ = format!("{h:?}");
    }

    #[test]
    fn test_getters() {
        let h = create_harmonic();
        h.path();
        h.typ();
        h.m();
        h.n();
        h.phase_average();
        h.len();

        assert_eq!(h.psip_data().ndim(), 1);
        assert_eq!(h.a_data().ndim(), 1);
        assert_eq!(h.phase_data().ndim(), 1);
    }

    #[test]
    fn test_data_extraction() {
        let h = create_harmonic();

        assert_eq!(h.psip_data().ndim(), 1);
        assert_eq!(h.a_data().ndim(), 1);
        assert_eq!(h.phase_data().ndim(), 1);
    }

    #[test]
    fn test_cache_update() {
        let h = create_harmonic();
        let mut acc = Accelerator::new();
        let mut cache = HarmonicCache::new();

        // dh_dt does not update the cache
        let (psip, theta, zeta) = (0.015, 0.0, 3.14);
        h.h(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dpsip(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dtheta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        h.dh_dzeta(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dt(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        assert_eq!(cache.misses, 1);
        assert_eq!(cache.hits, 3);
        let (psip, theta, zeta) = (0.01, 0.01, 3.15);
        h.dh_dpsip(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.h(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dtheta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        h.dh_dzeta(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dt(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        assert_eq!(cache.misses(), 2);
        assert_eq!(cache.hits(), 6);
    }

    #[test]
    fn test_spline_evals() {
        let h = create_harmonic();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        h.a(psip, &mut acc).unwrap();
        h.da_dpsip(psip, &mut acc).unwrap();
        h.phase(psip, &mut acc).unwrap();
    }
}
