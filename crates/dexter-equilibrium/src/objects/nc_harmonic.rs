//! Representation of a numerical equilibrium's single harmonic.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};

use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{EqError, EvalError};
use crate::{EquilibriumType, Harmonic, HarmonicCache};

/// Defines the calculation method of the phase `φ` in an [`NcHarmonic`].
#[derive(Default, Debug, Clone)]
pub enum PhaseMethod {
    /// Corresponds to `φ = 0`.
    #[default]
    Zero,
    /// Corresponds to `φ = const = the average of all the values of the phase array`.
    Average,
    /// Corresponds to `φ = const = the value of φ at the resonance m/n`.
    ///
    /// In the case that the resonance falls outside the last closed flux surface, or does not
    /// correspond to a valid q-factor value, it defaults to [`Zero`](PhaseMethod::Zero).
    Resonance,
    /// Interpolation over the phase array.
    Interpolation,
    /// Use a custom value for `φ = const`.
    Custom(f64),
}

impl PhaseMethod {
    /// The method to fall back to if a more complex method fails.
    pub(crate) fn fallback() -> Self {
        Self::Zero
    }
}

/// Used to create an [`NcHarmonic`].
#[non_exhaustive]
#[derive(Debug)]
pub struct NcHarmonicBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`] type (case-insensitive).
    interp_type: String,
    /// The `θ` frequency number.
    m: i64,
    /// The `θ` frequency number.
    n: i64,
    /// The calculation method of the phase `φ`.
    phase_method: PhaseMethod,
}

impl NcHarmonicBuilder {
    /// Creates a new [`NcHarmonicBuilder`] from a netCDF file at `path`, with spline of
    /// `interp1d_type` interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcHarmonicBuilder::new(&path, "steffen", 3, 2);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: &str, m: i64, n: i64) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type: interp_type.into(),
            m,
            n,
            phase_method: PhaseMethod::default(),
        }
    }

    /// Sets the phase `φ` calculation method.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcHarmonicBuilder::new(&path, "steffen", 3, 2)
    ///     .with_phase_method(PhaseMethod::Interpolation)
    ///     .build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    #[must_use]
    pub fn with_phase_method(mut self, method: PhaseMethod) -> Self {
        self.phase_method = method;
        self
    }

    /// Creates a new [`NcHarmonic`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// # let path = PathBuf::from("./netcdf.nc");
    /// let harmonic = NcHarmonicBuilder::new(&path, "cubic", 3, 2).build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcHarmonic`].
    pub fn build(self) -> Result<NcHarmonic, EqError> {
        NcHarmonic::build(self)
    }
}

/// Single perturbation harmonic from a netCDF file.
///
/// The harmonic has the form of `α(ψ/ψp) * cos(mθ-nζ+φ(ψ/ψp))`, where `α` and `φ` can be expressed
/// as functions of either or both `ψ`, `ψp`, and are calculated by interpolation over the
/// numerical data.
///
/// `φ` calculation can be further configured with the [`PhaseMethod`] helper struct.
///
/// Should be created with an [`NcHarmonicBuilder`].
#[non_exhaustive]
#[derive(Clone)]
pub struct NcHarmonic {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The interpolation type.
    interp_type: String,

    /// The harmonic's poloidal mode number `m`.
    m: i64,
    /// The harmonic's toroidal mode number `n`.
    n: i64,

    /// The harmonic as a function of `ψ`.
    psi_single: SingleNcHarmonic,
    /// The harmonic as a function of `ψp`.
    psip_single: SingleNcHarmonic,
}

impl NcHarmonic {
    /// Constructs an [`NcHarmonic`] from an [`NcHarmonicBuilder`].
    pub(crate) fn build(builder: NcHarmonicBuilder) -> Result<Self, EqError> {
        use crate::extract;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path.clone())?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);

        let psi_single = SingleNcHarmonic::build(&file, &builder, psi, "ψ".into())?;
        let psip_single = SingleNcHarmonic::build(&file, &builder, psip, "ψp".into())?;

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            m: builder.m,
            n: builder.n,
            psi_single,
            psip_single,
        })
    }

    /// Checks if the called evaluation method is defined, returning `Err()` if not.
    fn check_if_defined(state: &FluxCoordinateState, msg: &str) -> Result<(), EvalError> {
        if *state == FluxCoordinateState::Good {
            Ok(())
        } else {
            Err(EvalError::UndefinedEvaluation(msg.into()))
        }
    }
}

impl std::fmt::Debug for NcHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcHarmonic")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("interpolation type", &self.interp_type())
            .field("m", &self.m)
            .field("n", &self.n)
            .field("psi single", &self.psi_single)
            .field("psip single", &self.psip_single)
            .field("phase method", &self.psi_single.phase_method)
            .finish()
    }
}

/// Stores an [`NcHarmonic`]'s constant parameters and cached quantities.
#[derive(Debug, Clone)]
pub struct NcHarmonicCache {
    /// The number of cache hits.
    hits: usize,
    /// The number of cache misses.
    misses: usize,
    /// The harmonic's poloidal mode number `m`, casted to `f64` to use for evaluations.
    m: f64,
    /// The harmonic's toroidal mode number `n`, casted to `f64` to use for evaluations.
    n: f64,
    /// Ordered array with the intermediate values.
    ///
    /// cache = [
    ///     0 = flux
    ///     1 = theta
    ///     2 = zeta
    ///     3 = alpha
    ///     4 = dalpha
    ///     5 = phase
    ///     6 = modarg
    ///     7 = sin
    ///     8 = cos
    /// ].
    cache: [f64; 9],
    /// The Accelerator of the current flux coordinate.
    acc: Accelerator,
}

impl Default for NcHarmonicCache {
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            m: f64::NAN,
            n: f64::NAN,
            cache: [f64::NAN; 9],
            acc: Accelerator::new(),
        }
    }
}

impl HarmonicCache for NcHarmonicCache {
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) -> bool {
        #[expect(
            clippy::float_cmp,
            reason = "we want a cache hit only if all values are exactly equal"
        )]
        if (self.cache[0] == flux) && (self.cache[1] == theta) && (self.cache[2] == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) {
        self.cache[0] = flux;
        self.cache[1] = theta;
        self.cache[2] = zeta;

        self.cache[6] = (self.m * theta - self.n * zeta + self.cache[5]).rem_euclid(TAU);
        (self.cache[7], self.cache[8]) = self.cache[6].sin_cos();
    }

    harmonic_cache_counts_getter_impl!(NcHarmonicCache);
}

// Perform psi/psip debug assertions here since he have the extra information about the flux.
impl Harmonic for NcHarmonic {
    type Cache = NcHarmonicCache;

    fn psi_state(&self) -> FluxCoordinateState {
        self.psi_single.flux.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip_single.flux.state()
    }

    fn generate_cache(&self) -> Self::Cache {
        Self::Cache {
            m: self.m as f64,
            n: self.n as f64,
            cache: [f64::NAN; 9],
            ..Default::default()
        }
    }

    fn alpha_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "α(ψ)")?;
        self.psi_single.alpha(psi, theta, zeta, t, cache)
    }

    fn alpha_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "α(ψp)")?;
        self.psip_single.alpha(psip, theta, zeta, t, cache)
    }

    fn phase_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "φ(ψ)")?;
        self.psi_single.phase(psi, theta, zeta, t, cache)
    }

    fn phase_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "φ(ψp)")?;
        self.psip_single.phase(psip, theta, zeta, t, cache)
    }

    fn h_of_psi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "h(ψ)")?;
        self.psi_single.h(psi, theta, zeta, t, cache)
    }

    fn h_of_psip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "h(ψp)")?;
        self.psip_single.h(psip, theta, zeta, t, cache)
    }

    fn dh_dpsi(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dψ")?;
        self.psi_single.dh_dflux(psi, theta, zeta, t, cache)
    }

    fn dh_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dψp")?;
        self.psip_single.dh_dflux(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dtheta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dθ")?;
        self.psi_single.dh_dtheta(psi, theta, zeta, t, cache)
    }

    fn dh_of_psip_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dθ")?;
        self.psip_single.dh_dtheta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dzeta(
        &self,
        psi: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Self::check_if_defined(&self.psi_single.flux.state(), "dh(ψ)/dζ")?;
        self.psi_single.dh_dzeta(psi, theta, zeta, t, cache)
    }

    fn dh_of_psip_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Self::check_if_defined(&self.psip_single.flux.state(), "dh(ψp)/dζ")?;
        self.psip_single.dh_dzeta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dt(
        &self,
        psi: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(0.0)
    }

    fn dh_of_psip_dt(
        &self,
        psip: f64,
        _: f64,
        _: f64,
        _: f64,
        _: &mut Self::Cache,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(0.0)
    }
}

// Getters.
impl NcHarmonic {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);
    harmonic_mode_number_getter_impl!();

    /// Returns the `NcHarmonic`'s [`PhaseMethod`].
    ///
    /// Both flux coordinates have the same [`PhaseMethod`].
    #[must_use]
    pub fn phase_method(&self) -> PhaseMethod {
        self.psi_single.phase_method.clone()
    }

    /// Returns the average value of the phase arrays, if [`PhaseMethod`] is `Average`.
    #[must_use]
    pub fn phase_average(&self) -> Option<f64> {
        self.psi_single.phase_average
    }

    /// Returns the toroidal flux's value where the resonance is met, if [`PhaseMethod`] is
    /// `Resonance` and the resonance is in bounds.
    #[must_use]
    pub fn psi_phase_resonance(&self) -> Option<f64> {
        self.psi_single.phase_resonance
    }

    /// Returns the poloidal flux's value where the resonance is met, if [`PhaseMethod`] is
    /// `Resonance` and the resonance is in bounds.
    #[must_use]
    pub fn psip_phase_resonance(&self) -> Option<f64> {
        self.psip_single.phase_resonance
    }

    /// Returns the value of the last closed toroidal flux surface `ψ_last`.
    #[must_use]
    pub fn psi_last(&self) -> Option<f64> {
        self.psi_single.flux.last_value()
    }

    /// Returns the value of the last closed poloidal flux surface `ψp_last`.
    #[must_use]
    pub fn psip_last(&self) -> Option<f64> {
        self.psip_single.flux.last_value()
    }

    /// Returns the toroidal flux's values as a 1D array, if they exist.
    #[must_use]
    pub fn psi_array(&self) -> Option<Array1<f64>> {
        self.psi_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }
    /// Returns the poloidal flux's values as a 1D array, if they exist.
    #[must_use]
    pub fn psip_array(&self) -> Option<Array1<f64>> {
        self.psip_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }

    /// Returns the `α` values as a 1D array.
    #[must_use]
    pub fn alpha_array(&self) -> Array1<f64> {
        Array1::from_vec(self.psi_single.alpha_values.clone())
    }

    /// Returns the `φ` values as a 1D array.
    #[must_use]
    pub fn phase_array(&self) -> Array1<f64> {
        Array1::from_vec(self.psi_single.phase_values.clone())
    }
}

// ===============================================================================================
// ===============================================================================================

/// Representation of a numerical Harmonic, defined as a function of only one of the two flux
/// coordinates (and θ, ζ, t).
///
/// Splitting an [`NcHarmonic`] into two `SingleNcHarmonics` makes the code much clearer.
///
/// Both harmonics end up identical, with the only difference being the corresponding [`NcFlux`].
#[non_exhaustive]
struct SingleNcHarmonic {
    /// The interpolation type.
    interp_type: String,
    /// The current flux coordinate.
    flux: NcFlux,
    /// `ψ` or `ψp`, to be used in [`EqError::UndefinedEvaluation`] message.
    which: Box<str>,
    /// The harmonic's poloidal mode number `m`.
    m: i64,
    /// The harmonic's toroidal mode number `n`.
    n: i64,
    /// The phase calculation method.
    phase_method: PhaseMethod,
    /// The phase values' average, if `phase_method` is `Average`.
    phase_average: Option<f64>,
    /// The phase's value at the resonance, if `phase_method` is `Resonance`.
    phase_resonance: Option<f64>,

    /// The amplitude values.
    alpha_values: Vec<f64>,
    /// The phase values.
    phase_values: Vec<f64>,
    /// The `α(flux)` interpolator.
    alpha_interp: Option<DynInterpolation<f64>>,
    /// The `φ(flux)` interpolator.
    phase_interp: Option<DynInterpolation<f64>>,
}

// Unforturately we must rebuild the interpolators, since they are trait objects.
impl Clone for SingleNcHarmonic {
    fn clone(&self) -> Self {
        Self {
            interp_type: self.interp_type.clone(),
            flux: self.flux.clone(),
            which: self.which.clone(),
            m: self.m,
            n: self.n,
            phase_method: self.phase_method.clone(),
            phase_average: self.phase_average,
            phase_resonance: self.phase_resonance,
            alpha_values: self.alpha_values.clone(),
            phase_values: self.phase_values.clone(),
            alpha_interp: make_interp_type(&self.interp_type)
                .expect("Already built once, cannot fail")
                .build(self.flux.uvalues(), &self.alpha_values)
                .ok(),
            phase_interp: make_interp_type(&self.interp_type)
                .expect("Already built once, cannot fail")
                .build(self.flux.uvalues(), &self.phase_values)
                .ok(),
        }
    }
}

// Creation.
impl SingleNcHarmonic {
    /// Creates a `SingleNcHarmonic` from a [`NcFlux`].
    ///
    /// The object always exists, even if the corresponding flux coordinate does not exist, so all
    /// logic and error handling is handle here.
    pub(crate) fn build(
        file: &netcdf::File,
        builder: &NcHarmonicBuilder,
        flux: NcFlux,
        which: Box<str>,
    ) -> Result<Self, EqError> {
        use crate::extract;

        let (alpha_data, phase_data) = extract::harmonic_arrays(file, builder.m, builder.n)?;
        let alphas = alpha_data.to_vec();
        let phases = phase_data.to_vec();

        debug_assert_all_finite_values(&alphas);
        debug_assert_all_finite_values(&phases);

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let alpha_interp = match flux.state() {
            Good => Some(make_interp_type(&builder.interp_type)?.build(flux.uvalues(), &alphas)?),
            _ => None,
        };
        #[rustfmt::skip]
        let phase_interp = match flux.state() {
            Good => Some(make_interp_type(&builder.interp_type)?.build(flux.uvalues(), &phases)?),
            _ => None,
        };

        let _m = builder.m as f64;
        let _n = builder.n as f64;

        let uninit_self = Self {
            interp_type: builder.interp_type.clone(),
            flux,
            which,
            m: builder.m,
            n: builder.n,
            phase_method: builder.phase_method.clone(),
            phase_average: None,
            phase_resonance: None,
            alpha_values: alphas,
            phase_values: phases,
            alpha_interp,
            phase_interp,
        };
        let init_self = uninit_self.resolve_phase_method(file);
        Ok(init_self)
    }

    /// Calculates the phase method. Should only be used on initialization.
    fn resolve_phase_method(self, file: &netcdf::File) -> Self {
        let phase_method: PhaseMethod;
        let mut phase_average: Option<f64> = None;
        let mut phase_resonance: Option<f64> = None;

        match self.phase_method {
            PhaseMethod::Zero => phase_method = PhaseMethod::Zero,
            PhaseMethod::Interpolation => phase_method = PhaseMethod::Interpolation,
            PhaseMethod::Custom(phase) => phase_method = PhaseMethod::Custom(phase),
            PhaseMethod::Average => {
                phase_method = PhaseMethod::Average;
                phase_average = Array1::from_vec(self.phase_values.clone()).mean();
            }
            PhaseMethod::Resonance => match self.find_resonance_phase(file) {
                Some(value) => {
                    phase_method = PhaseMethod::Resonance;
                    phase_resonance = Some(value);
                }
                None => {
                    phase_method = PhaseMethod::fallback();
                }
            },
        }

        Self {
            phase_method,
            phase_average,
            phase_resonance,
            ..self
        }
    }

    /// Calculates the phase's value at the resonance.
    fn find_resonance_phase(&self, _: &netcdf::File) -> Option<f64> {
        unimplemented!()
    }
}

/// External cache update.
impl SingleNcHarmonic {
    /// Updates the interpolated values, since they cannot take place inside the cache without
    /// messing up the whole structure.
    ///
    /// Should always be called together with `cache.update()`, and `update_cache_interps()` should
    /// always be called firsts, since `cache` uses the phase value.
    #[rustfmt::skip]
    fn update_cache_interps(
        &self,
        flux: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<(), EvalError> {
        let acc = &mut cache.acc;

        cache.cache[3] = match self.alpha_interp.as_ref() {
            Some(interp) => interp.eval(self.flux.uvalues(), &self.alpha_values, flux, acc)?,
            None => return Err(EvalError::UndefinedEvaluation("α(flux)".into())),
        };
        cache.cache[4] = match self.alpha_interp.as_ref() {
            Some(interp) => interp.eval_deriv(self.flux.uvalues(), &self.alpha_values, flux, acc)?,
            None => return Err(EvalError::UndefinedEvaluation("da(flux)/dflux".into())),
        };
        cache.cache[5] = match self.phase_method {
            PhaseMethod::Zero => 0.0,
            PhaseMethod::Average => self.phase_average.expect("Exists"),
            PhaseMethod::Resonance => self.phase_resonance.expect("Exists"),
            PhaseMethod::Custom(custom_phase) => custom_phase,
            PhaseMethod::Interpolation => {
                match self.phase_interp.as_ref() {
                    Some(interp) => interp.eval(self.flux.uvalues(), &self.phase_values, flux, acc)?,
                    None => return Err(EvalError::UndefinedEvaluation("φ(flux)".into())),
                }
            }
        };
        Ok(())
    }
}

/// Intermediate Interpolations. This is effectively where the [`Harmonic`] trait is implemented.
impl SingleNcHarmonic {
    /// Calculates the single harmonic's amplitude.
    fn alpha(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache[3]))
    }

    /// Calculates the single harmonic's phase.
    fn phase(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache[5]))
    }

    /// Calculates the single harmonic's value.
    fn h(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache[3] * cache.cache[8]))
    }

    /// Calculates the single harmonic's derivative with respect to the current flux coordinate.
    fn dh_dflux(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(cache.cache[4] * cache.cache[8]))
    }

    /// Calculates the single harmonic's derivative with respect to theta.
    fn dh_dtheta(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(
            -cache.m * cache.cache[3] * cache.cache[7]
        ))
    }

    /// Calculates the single harmonic's derivative with respect to zeta.
    fn dh_dzeta(
        &self,
        flux: f64,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut NcHarmonicCache,
    ) -> Result<f64, EvalError> {
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t);
        };
        Ok(debug_assert_is_finite!(
            cache.n * cache.cache[3] * cache.cache[7]
        ))
    }
}

impl std::fmt::Debug for SingleNcHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SingleNcHarmonic")
            .field("NcFlux", &self.flux)
            .field("phase method", &self.phase_method)
            .field("phase average", &self.phase_average)
            .finish()
    }
}

// ===============================================================================================

#[cfg(test)]
mod test_utils {
    use super::*;

    pub(super) fn create_nc_harmonic_builder(path: &str) -> NcHarmonicBuilder {
        let path = PathBuf::from(path);
        NcHarmonicBuilder::new(&path, "steffen", 3, 2)
    }

    pub(super) fn create_nc_harmonic(path: &str) -> NcHarmonic {
        create_nc_harmonic_builder(path)
            .with_phase_method(PhaseMethod::Zero)
            .build()
            .unwrap()
    }
}

#[cfg(test)]
mod phase_methods {
    use crate::extract::{POLOIDAL_TEST_NETCDF_PATH, TEST_NETCDF_PATH};
    use approx::assert_relative_eq;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn zero_phase_method() {
        use PhaseMethod::Zero;
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            .with_phase_method(Zero)
            .build()
            .unwrap();
        let c = &mut h.generate_cache();

        assert!(matches!(h.phase_method(), Zero));
        assert!(h.psi_single.phase_average.is_none());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_eq!(h.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), 0.0);
        assert_eq!(h.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), 0.0);
    }

    #[test]
    fn average_phase_method() {
        use PhaseMethod::Average;
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            .with_phase_method(Average)
            .build()
            .unwrap();
        let c = &mut h.generate_cache();

        let expected = h.phase_array().mean().unwrap();

        assert!(matches!(h.phase_method(), Average));
        assert!(h.psi_single.phase_average.is_some());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_eq!(h.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), expected);
        assert_eq!(h.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), expected);
    }

    #[test]
    #[ignore = "re-write"]
    fn resonance_phase_method() {
        use PhaseMethod::Resonance;
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let h = NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .with_phase_method(Resonance)
            .build()
            .unwrap();
        let c = &mut h.generate_cache();

        assert!(h.psip_single.phase_average.is_none());

        // m=3 and n=2, we expect a resonance at q=n/m=2/3=0.666, which is inbounds in our
        // poloidal_test_netcdf. However, only ψp is a good coordinate.

        // From the create_test_netcdf script:
        let q_res = h.n() as f64 / h.m() as f64;
        let psip_res = TAU * (TAU * q_res).cos(); // we defined ψ = sin(2πψp)
        let expected = h.phase_of_psip(psip_res, 0.1, 0.1, 0.1, c).unwrap();

        assert!(matches!(h.phase_method(), Resonance));
        assert!(h.psip_single.phase_average.is_none());
        assert_relative_eq!(
            h.psip_single.phase_resonance.unwrap(),
            expected,
            epsilon = 1e-10
        );
    }

    #[test]
    fn interpolation_phase_method() {
        use PhaseMethod::Interpolation;
        let h = create_nc_harmonic_builder(POLOIDAL_TEST_NETCDF_PATH)
            .with_phase_method(Interpolation)
            .build()
            .unwrap();
        let c = &mut h.generate_cache();

        // Calculated with a stable version, on the same dataset
        let expected = 0.8414720460746888;

        assert!(matches!(h.phase_method(), Interpolation));
        assert!(h.psi_single.phase_average.is_none());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_relative_eq!(
            h.phase_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap(),
            expected,
            epsilon = 1e-5 // unsure why there is this small difference here.
        );
    }

    #[test]
    fn custom_phase_method() {
        use PhaseMethod::Custom;
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            .with_phase_method(Custom(10.0))
            .build()
            .unwrap();
        let c = &mut h.generate_cache();

        assert!(matches!(h.phase_method(), Custom(10.0)));
        assert!(h.psi_single.phase_average.is_none());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_eq!(h.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
        assert_eq!(h.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
    }

    #[test]
    #[ignore = "re-write"]
    fn fallback_phase_method() {
        let path = PathBuf::from(TEST_NETCDF_PATH);
        let h = NcHarmonicBuilder::new(&path, "steffen", 2, 2)
            .with_phase_method(PhaseMethod::Resonance)
            .build()
            .unwrap();
        dbg!(&h);

        // q=1 on our test_netcdf, so no resonance should be found for m=3, n=2, since q_res = 2/3

        assert!(h.psip_single.phase_average.is_none());
        assert!(h.psip_single.phase_resonance.is_none());
        assert!(matches!(h.phase_method(), PhaseMethod::Zero));
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let harmonic = create_nc_harmonic(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(harmonic.psi_state(), FluxCoordinateState::Good);
        assert_eq!(harmonic.psip_state(), FluxCoordinateState::Bad);
        assert!(harmonic.psi_single.alpha_interp.is_some());
        assert!(harmonic.psi_single.phase_interp.is_some());
        assert!(harmonic.psip_single.alpha_interp.is_none());
        assert!(harmonic.psip_single.phase_interp.is_none());

        assert!(harmonic.psi_array().is_some());
        assert!(harmonic.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psi_evals() {
        let harmonic = create_nc_harmonic(TOROIDAL_TEST_NETCDF_PATH);
        let c = &mut harmonic.generate_cache();
        assert!(harmonic.alpha_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.phase_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.h_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_dpsi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psi_dtheta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psi_dzeta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psi_dt(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
    }

    #[test]
    #[rustfmt::skip]
    fn bad_psip_evals() {
        let h = create_nc_harmonic(TOROIDAL_TEST_NETCDF_PATH);
        let c = &mut h.generate_cache();

        use EvalError::UndefinedEvaluation as err;
        assert!(matches!(h.alpha_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.phase_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.h_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_dpsip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psip_dtheta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psip_dzeta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psip_dt(0.1, 0.1, 0.1, 0.1, c), Ok(0.0)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let harmonic = create_nc_harmonic(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(harmonic.psi_state(), FluxCoordinateState::Bad);
        assert_eq!(harmonic.psip_state(), FluxCoordinateState::Good);
        assert!(harmonic.psi_single.alpha_interp.is_none());
        assert!(harmonic.psi_single.phase_interp.is_none());
        assert!(harmonic.psip_single.alpha_interp.is_some());
        assert!(harmonic.psip_single.phase_interp.is_some());

        assert!(harmonic.psi_array().is_some());
        assert!(harmonic.psip_array().is_some());
    }

    #[test]
    #[rustfmt::skip]
    fn good_psip_evals() {
        let harmonic = create_nc_harmonic(POLOIDAL_TEST_NETCDF_PATH);
        let c = &mut harmonic.generate_cache();
        assert!(harmonic.alpha_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.phase_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.h_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_dpsip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psip_dtheta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psip_dzeta(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
        assert!(harmonic.dh_of_psip_dt(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
    }

    #[test]
    #[rustfmt::skip]
    fn bad_psi_evals() {
        let h = create_nc_harmonic(POLOIDAL_TEST_NETCDF_PATH);
        let c = &mut h.generate_cache();

        use EvalError::UndefinedEvaluation as err;
        assert!(matches!(h.alpha_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.phase_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.h_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_dpsi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psi_dtheta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psi_dzeta(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.dh_of_psi_dt(0.1, 0.1, 0.1, 0.1, c), Ok(0.0)));
    }
}

#[cfg(test)]
mod nc_harmonic_cache {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    #[allow(unused_results)]
    fn counts() {
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            // other methods don't touch the cache at all
            .with_phase_method(PhaseMethod::Interpolation)
            .build()
            .unwrap();
        let c = &mut h.generate_cache();
        let t = 0.0; // not checked

        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 0);

        h.phase_of_psi(0.01, 0.1, 0.1, t, c).unwrap();
        assert_eq!(c.hits(), 0);
        assert_eq!(c.misses(), 1);

        h.dh_of_psi_dtheta(0.01, 0.1, 0.1, t, c).unwrap();
        assert_eq!(c.hits(), 1);
        assert_eq!(c.misses(), 1);

        let psi = 0.1;
        let theta = 3.14;
        let zeta = 1.0;
        h.alpha_of_psi(psi, theta, zeta, t, c).unwrap();
        h.phase_of_psi(psi, theta, zeta, t, c).unwrap();
        h.h_of_psi(psi, theta, zeta, t, c).unwrap();
        h.dh_of_psi_dtheta(psi, theta, zeta, t, c).unwrap();
        h.dh_of_psi_dzeta(psi, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 5);
        assert_eq!(c.misses(), 2);

        h.dh_dpsi(psi / 2.0, theta, zeta, t, c).unwrap();

        assert_eq!(c.hits(), 5);
        assert_eq!(c.misses(), 3);

        // dbg!(&c); // Accelerator counts should also match
    }
}
