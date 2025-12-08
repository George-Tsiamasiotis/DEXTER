//! Representation of an perturbation 's single harmonic.

use std::path::PathBuf;

use common::array1D_getter_impl;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};

use ndarray::Array1;

use crate::Result;
use crate::{Flux, Length, Radians};
use crate::{Harmonic, HarmonicCache};

/// Defines the calculation method of the phase `φ(ψp)`.
#[derive(Default, Debug, Clone)]
pub enum PhaseMethod {
    /// φ(ψp) = 0.
    Zero,
    /// φ = const = the average of all the values of the phase array.
    Average,
    /// φ = const = the value of φ(ψp) at the resonance `m/n`.
    #[default]
    Resonance,
    /// Interpolate over the phase array.
    Interpolation,
    /// Use a custom value for φ = const.
    Custom(f64),
}

/// Used to create an [`NcHarmonic`].
pub struct NcHarmonicBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,
    /// The `θ` frequency number.
    m: i64,
    /// The `θ` frequency number.
    n: i64,
    /// The calculation method of the phase `φ(ψp)`.
    phase_method: PhaseMethod,
}

impl NcHarmonicBuilder {
    /// Creates a new [`NcHarmonicBuilder`] from a netCDF file at `path`, with spline of `typ`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, "cubic");
    /// ```
    pub fn new(path: &PathBuf, typ: &str, m: i64, n: i64) -> Self {
        Self {
            path: path.clone(),
            typ: typ.into(),
            m,
            n,
            phase_method: PhaseMethod::default(),
        }
    }

    /// Sets the phase `φ(ψp)` calculation method.
    pub fn phase_method(mut self, method: PhaseMethod) -> Self {
        self.phase_method = method;
        self
    }

    /// Creates a new [`NcHarmonic`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use equilibrium::*;
    /// # use std::path::PathBuf;
    /// #
    /// # fn main() -> Result<()> {
    /// let path = PathBuf::from("../data/stub_netcdf.nc");
    /// let harmonic = NcHarmonicBuilder::new(&path, "cubic").build()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build(self) -> Result<NcHarmonic> {
        NcHarmonic::build(self)
    }
}

// ===============================================================================================

/// Single perturbation harmonic reconstructed from a netCDF file.
///
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ(ψp))`, where `α(ψp)` is calculated by
/// interpolation over numerical data, and `φ(ψp)` is calculated as defined by [`PhaseMethod`].
pub struct NcHarmonic {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`Interpolation type`], in case-insensitive string format.
    ///
    /// [`Interpolation type`]: ../rsl_interpolation/trait.InterpType.html#implementors
    typ: String,

    /// The `θ` frequency number, cast to f64 to be used in calculations.
    pub(crate) _m: f64,
    /// The `ζ` frequency number, cast to f64 to be used in calculations.
    pub(crate) _n: f64,
    /// The calculation method of the phase `φ(ψp)`.
    pub(crate) phase_method: PhaseMethod,
    /// The average value of the phase array.
    pub(crate) phase_average: Option<f64>,
    /// The value of the phase at the resonance ψp = m/n.
    pub(crate) phase_resonance: Option<f64>,

    /// The `ψp` data array.
    psip_data: Vec<Flux>,
    /// The `α` data array.
    alpha_data: Vec<Length>,
    /// The `φ` data array.
    phase_data: Vec<Radians>,

    /// Interpolator over the `α` values, as a function of ψp.
    alpha_interp: DynInterpolation<f64>,
    /// Interpolator over the `α` values, as a function of ψp.
    phase_interp: DynInterpolation<f64>,
}

/// Creation
impl NcHarmonic {
    /// Constructs an [`NcHarmonic`] from [`NcHarmonicBuilder`].
    pub(crate) fn build(builder: NcHarmonicBuilder) -> Result<Self> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let f = open(&path)?;

        let psip_data = extract_1d_array(&f, PSIP_NORM)?
            .as_standard_layout()
            .to_vec();
        let (alpha_data, phase_data) = extract_harmonic_arrays(&f, builder.m, builder.n)?;
        let alpha_data = alpha_data.to_vec();
        let phase_data = phase_data.to_vec();

        let alpha_interp = make_interp_type(&builder.typ)?.build(&psip_data, &alpha_data)?;
        let phase_interp = make_interp_type(&builder.typ)?.build(&psip_data, &phase_data)?;

        let phase_average = match builder.phase_method {
            PhaseMethod::Average => Some(
                // If `phase_data` was empty, `extract_1d_array` would have failed.
                Array1::from(phase_data.clone())
                    .mean()
                    .expect("array is non-empty"),
            ),
            _ => None,
        };

        let phase_resonance = match builder.phase_method {
            PhaseMethod::Resonance => {
                let mut acc = Accelerator::new();
                Some(phase_interp.eval(
                    &psip_data,
                    &phase_data,
                    (builder.m as f64) / (builder.n as f64),
                    &mut acc,
                )?)
            }
            _ => None,
        };

        Ok(Self {
            path: path.to_owned(),
            typ: builder.typ,
            _m: builder.m as f64,
            _n: builder.n as f64,
            phase_method: builder.phase_method,
            psip_data,
            alpha_data,
            phase_average,
            phase_resonance,
            phase_data,
            alpha_interp,
            phase_interp,
        })
    }
}

impl Harmonic for NcHarmonic {
    fn h(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha() * cache.cos())
    }

    fn dh_dpsip(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.dalpha() * cache.cos())
    }

    fn dh_dtheta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha() * (-self._m) * cache.sin())
    }

    fn dh_dzeta(
        &self,
        psip: Flux,
        theta: Radians,
        zeta: Radians,
        cache: &mut Box<dyn HarmonicCache>,
        acc: &mut Accelerator,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha() * self._n * cache.sin())
    }

    fn a(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .alpha_interp
            .eval(&self.psip_data, &self.alpha_data, psip, acc)?)
    }

    fn da_dpsip(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .alpha_interp
            .eval_deriv(&self.psip_data, &self.alpha_data, psip, acc)?)
    }

    fn phase(&self, psip: Flux, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .phase_interp
            .eval(&self.psip_data, &self.alpha_data, psip, acc)?)
    }
}

/// Getters
impl NcHarmonic {
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
        self.psip_data.len()
    }

    /// Returns the poloidal mode number `m`.
    pub fn m(&self) -> i64 {
        self._m as i64
    }

    /// Returns the poloidal mode number `m`.
    pub fn n(&self) -> i64 {
        self._n as i64
    }

    /// Returns the [`NcHarmonic`]'s phase calculation method.
    pub fn phase_method(&self) -> PhaseMethod {
        self.phase_method.clone()
    }

    /// Returns the [`NcHarmonic`]'s phase average.
    ///
    /// Returns `None` if the [`NcHarmonic`]'s [`PhaseMethod`] is not [`PhaseMethod::Average`].
    pub fn phase_average(&self) -> Option<Radians> {
        self.phase_average
    }

    /// Returns the [`NcHarmonic`]'s phase value at the resonance φ(ψp) = m/n.
    ///
    /// Returns `None` if the [`NcHarmonic`]'s [`PhaseMethod`] is not [`PhaseMethod::Resonance`].
    pub fn phase_resonance(&self) -> Option<Radians> {
        self.phase_resonance
    }

    array1D_getter_impl!(psip_data, psip_data, Flux);
    array1D_getter_impl!(a_data, alpha_data, Length);
    array1D_getter_impl!(phase_data, phase_data, Radians);
}

impl Clone for NcHarmonic {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            typ: self.typ.clone(),
            _m: self._m,
            _n: self._n,
            phase_method: self.phase_method.clone(),
            phase_average: self.phase_average,
            phase_resonance: self.phase_resonance,
            psip_data: self.psip_data.clone(),
            alpha_data: self.alpha_data.clone(),
            phase_data: self.phase_data.clone(),
            alpha_interp: make_interp_type(&self.typ)
                .expect("already created with the same typ")
                .build(&self.psip_data, &self.alpha_data)
                .expect("already created with the same arrays"),
            phase_interp: make_interp_type(&self.typ)
                .expect("already created with the same typ")
                .build(&self.psip_data, &self.phase_data)
                .expect("already created with the same arrays"),
        }
    }
}

impl std::fmt::Debug for NcHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcHarmonic")
            .field("path", &self.path())
            .field("typ", &self.typ())
            .field("m", &self.m())
            .field("n", &self.n())
            .field("phase_method", &self.phase_method)
            .field("phase_average", &self.phase_average)
            .field("phase_resonance", &self.phase_resonance)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::extract::STUB_NETCDF_PATH;
    use crate::*;

    fn create_nc_harmonic() -> NcHarmonic {
        let path = PathBuf::from(STUB_NETCDF_PATH);
        let typ = "steffen";
        NcHarmonicBuilder::new(&path, typ, 3, 2).build().unwrap()
    }

    #[test]
    fn test_creation() {
        let h = create_nc_harmonic();
        let _ = h.clone();
        let _ = format!("{h:?}");
    }

    #[test]
    fn test_getters() {
        let h = create_nc_harmonic();
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
        let h = create_nc_harmonic();

        assert_eq!(h.psip_data().ndim(), 1);
        assert_eq!(h.a_data().ndim(), 1);
        assert_eq!(h.phase_data().ndim(), 1);
    }

    #[test]
    fn test_cache_update() {
        let h = create_nc_harmonic();
        let mut acc = Accelerator::new();
        let mut cache: Box<dyn HarmonicCache> = Box::new(NcHarmonicCache::new());

        // dh_dt does not update the cache
        let (psip, theta, zeta) = (0.015, 0.0, 3.14);
        h.h(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dpsip(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dtheta(psip, theta, zeta, &mut cache, &mut acc)
            .unwrap();
        h.dh_dzeta(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        h.dh_dt(psip, theta, zeta, &mut cache, &mut acc).unwrap();
        assert_eq!(cache.misses(), 1);
        assert_eq!(cache.hits(), 3);
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
        let h = create_nc_harmonic();
        let mut acc = Accelerator::new();

        let psip = 0.015;
        h.a(psip, &mut acc).unwrap();
        h.da_dpsip(psip, &mut acc).unwrap();
        h.phase(psip, &mut acc).unwrap();
    }
}
