//! Representation of a numerical equilibrium's single harmonic.

use crate::eval::HarmonicCache;
use crate::flux::{NcFlux, NcFluxState};
use crate::{EqError, EquilibriumType, Harmonic, Result};
use crate::{
    equilibrium_type_getter_impl, harmonic_cache_counts_getter_impl,
    harmonic_mode_number_getter_impl, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};

use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

/// Defines the calculation method of the phase `φ` in a Numerical [`Harmonic`].
#[derive(Default, Debug, Clone)]
pub enum PhaseMethod {
    /// Corresponds to `φ = 0`.
    Zero,
    /// Corresponds to `φ = const = the average of all the values of the phase array`.
    Average,
    /// Corresponds to `φ = const = the value of φ at the resonance m/n`.
    ///
    /// In the case that the resonance falls outside the wall, or does not correspond to a valid
    /// q-factor value, it defaults to [`Zero`](PhaseMethod::Zero).
    #[default]
    Resonance,
    /// Interpolation over the phase array.
    Interpolation,
    /// Use a custom value for φ = const.
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
    pub fn build(self) -> Result<NcHarmonic> {
        NcHarmonic::build(self)
    }
}

/// Single perturbation harmonic from a netCDF file.
///
/// The harmonic has the form of `α(ψ/ψp) * cos(mθ-nζ+φ(ψ/ψp))`, where `α` and `φ` can be expressed
/// as functions of either or both `ψ`, `ψp`, and are calculated by interpolation over the
/// numerical data,
///
/// `φ` calculation can be further configured with the [`PhaseMethod`] helper struct.
///
/// Should be created with an [`NcHarmonicBuilder`].
#[non_exhaustive]
#[derive(Clone)]
pub struct NcHarmonic {
    path: PathBuf,
    netcdf_version: semver::Version,

    equilibrium_type: EquilibriumType,
    interp_type: String,

    m: i64,
    n: i64,

    psi_single: SingleNcHarmonic,
    psip_single: SingleNcHarmonic,
}

// Creation
impl NcHarmonic {
    pub(crate) fn build(builder: NcHarmonicBuilder) -> Result<Self> {
        use crate::extract::*;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path.clone())?;
        let f = open(&path)?;
        let netcdf_version = extract_version(&f)?;

        let psi = NcFlux::toroidal(&f);
        let psip = NcFlux::poloidal(&f);

        let psi_single = SingleNcHarmonic::build(&f, &builder, psi, "ψ".into())?;
        let psip_single = SingleNcHarmonic::build(&f, &builder, psip, "ψp".into())?;

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
}

#[rustfmt::skip]
impl Harmonic for NcHarmonic
{
    type Cache = NcHarmonicCache;

    fn ampl_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.ampl(psi, theta, zeta, t, cache)
    }

    fn ampl_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.ampl(psip, theta, zeta, t, cache)
    }

    fn phase_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.phase(psi, theta, zeta, t, cache)
    }

    fn phase_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.phase(psip, theta, zeta, t, cache)
    }

    fn h_of_psi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.h(psi, theta, zeta, t, cache)
    }

    fn h_of_psip(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.h(psip, theta, zeta, t, cache)
    }

    fn dh_dpsi(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.dh_dflux(psi, theta, zeta, t, cache)
    }

    fn dh_dpsip(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.dh_dflux(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dtheta(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.dh_dtheta(psi, theta, zeta, t, cache)
    }

    fn dh_of_psip_dtheta(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.dh_dtheta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dzeta(&self, psi: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psi_single.dh_dzeta(psi, theta, zeta, t, cache)
    }

    fn dh_of_psip_dzeta(&self, psip: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.psip_single.dh_dzeta(psip, theta, zeta, t, cache)
    }

    fn dh_of_psi_dt(&self, _: f64, _: f64, _: f64, _: f64, _: &mut NcHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }

    fn dh_of_psip_dt(&self, _: f64, _: f64, _: f64, _: f64, _: &mut NcHarmonicCache) -> Result<f64> {
        Ok(0.0)
    }
}

// Getters
impl NcHarmonic {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);
    harmonic_mode_number_getter_impl!();

    /// Returns the NcHarmonic's [`PhaseMethod`].
    ///
    /// Both flux coordinates have the same [`PhaseMethod`].
    pub fn phase_method(&self) -> PhaseMethod {
        self.psi_single.phase_method.clone()
    }

    /// Returns the average value of the phase arrays, if [`PhaseMethod`] is `Average`.
    pub fn phase_average(&self) -> Option<f64> {
        self.psi_single.phase_average
    }

    /// Returns the toroidal flux's value where the resonance is met, if [`PhaseMethod`] is
    /// `Resonance` and the resonance is in bounds.
    pub fn psi_phase_resonance(&self) -> Option<f64> {
        self.psi_single.phase_resonance
    }

    /// Returns the poloidal flux's value where the resonance is met, if [`PhaseMethod`] is
    /// `Resonance` and the resonance is in bounds.
    pub fn psip_phase_resonance(&self) -> Option<f64> {
        self.psip_single.phase_resonance
    }

    /// Returns the toroidal flux's value at the wall `ψ_wall`.
    pub fn psi_wall(&self) -> Option<f64> {
        self.psi_single.flux.wall_value()
    }

    /// Returns the poloidal flux's value at the wall `ψp_wall`.
    pub fn psip_wall(&self) -> Option<f64> {
        self.psip_single.flux.wall_value()
    }

    /// Returns the toroidal flux's state.
    pub fn psi_state(&self) -> NcFluxState {
        self.psi_single.flux.state.clone()
    }
    /// Returns the poloidal flux's state.
    pub fn psip_state(&self) -> NcFluxState {
        self.psip_single.flux.state.clone()
    }

    /// Returns the toroidal flux's values as a 1D array, if they exist.
    pub fn psi_array(&self) -> Option<Array1<f64>> {
        self.psi_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }
    /// Returns the poloidal flux's values as a 1D array, if they exist.
    pub fn psip_array(&self) -> Option<Array1<f64>> {
        self.psip_single
            .flux
            .values()
            .map(|values| Array1::from(Vec::from(values)))
    }

    /// Returns the `α` values as a 1D array.
    pub fn alpha_array(&self) -> Array1<f64> {
        Array1::from_vec(self.psi_single.alpha_values.clone())
    }

    /// Returns the `φ` values as a 1D array.
    pub fn phase_array(&self) -> Array1<f64> {
        Array1::from_vec(self.psi_single.phase_values.clone())
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

// ===============================================================================================
// ===============================================================================================

/// Representation of a numerical Harmonic, defined as a function of only one of the two flux
/// coordinates (and θ, ζ, t).
///
/// Splitting an [`NcHarmonic`] into two SingleFluxNcHarmonics makes the code much clearer.
///
/// Both harmonics end up identical, with the only difference being the corresponding [`NcFlux`].
#[non_exhaustive]
pub(crate) struct SingleNcHarmonic {
    interp_type: String,
    flux: NcFlux,
    /// `ψ` or `ψp`, to be used in [`EqError::UndefinedEvaluation`] message.
    which: Box<str>,
    m: i64,
    n: i64,
    phase_method: PhaseMethod,
    phase_average: Option<f64>,
    phase_resonance: Option<f64>,

    alpha_values: Vec<f64>,
    phase_values: Vec<f64>,
    alpha_interp: Option<DynInterpolation<f64>>,
    phase_interp: Option<DynInterpolation<f64>>,

    _m: f64,
    _n: f64,
    // [m, n]
    _cache_args: [f64; 2],
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
            alpha_interp: Some(
                make_interp_type(&self.interp_type)
                    .unwrap()
                    .build(self.flux.uvalues(), &self.alpha_values)
                    .unwrap(),
            ),
            phase_interp: Some(
                make_interp_type(&self.interp_type)
                    .unwrap()
                    .build(self.flux.uvalues(), &self.phase_values)
                    .unwrap(),
            ),
            _m: self._m,
            _n: self._n,
            _cache_args: self._cache_args,
        }
    }
}

// Creation
impl SingleNcHarmonic {
    /// Creates a SingleNcHarmonic from a [`NcFlux`].
    ///
    /// The object always exists, even if the corresponding flux coordinate does not exist, so all
    /// logic and error handling is handle here.
    pub(crate) fn build(
        f: &netcdf::File,
        builder: &NcHarmonicBuilder,
        flux: NcFlux,
        which: Box<str>,
    ) -> Result<Self> {
        use crate::extract::*;

        let (alpha_data, phase_data) = extract_harmonic_arrays(f, builder.m, builder.n)?;
        let alpha_values = alpha_data.to_vec();
        let phase_values = phase_data.to_vec();

        // Create interpolators, if possible
        use NcFluxState::Good;
        #[rustfmt::skip]
        let alpha_interp = match flux.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(flux.uvalues(), &alpha_values)?),
            _ => None,
        };
        #[rustfmt::skip]
        let phase_interp = match flux.state {
            Good => Some(make_interp_type(&builder.interp_type)?.build(flux.uvalues(), &phase_values)?),
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
            alpha_values,
            phase_values,
            alpha_interp,
            phase_interp,
            _m,
            _n,
            _cache_args: [_m, _n],
        };
        let init_self = uninit_self.resolve_phase_method(f);
        Ok(init_self)
    }

    /// Calculates the phase method. Should only be used on initialization.
    fn resolve_phase_method(self, f: &netcdf::File) -> Self {
        let phase_method: PhaseMethod;
        let mut phase_average: Option<f64> = None;
        let mut phase_resonance: Option<f64> = None;

        use PhaseMethod::*;
        match self.phase_method {
            Zero => phase_method = Zero,
            Interpolation => phase_method = Interpolation,
            Custom(phase) => phase_method = Custom(phase),
            Average => {
                phase_method = Average;
                phase_average = Array1::from_vec(self.phase_values.clone()).mean();
            }
            Resonance => match self.find_resonance_phase(f) {
                Some(value) => {
                    phase_method = Resonance;
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

    /// Attempts to find the value of the phase `φ` at the resonance. If it fails, it fallback to
    /// [`PhaseMethod::fallback`].
    ///
    /// For the attempt to succeed, it is necessary that:
    ///     - the q-factor values exist
    ///     - the corresponding flux coordinate's values exist
    ///     - the q-factor values are monotonic, so there is no ambiguity on the q->flux
    ///         interpolation.
    ///     - the q-factor of the resonance is outside the bounds of the q-factor values
    ///
    /// The 3rd assumption can be eliminated by assuming the q-factor values and phase values
    /// have the same length, and just picking the phase value with the same index as the q-factor
    /// element closer to the resonance. This introduces a small error, but might be the better way
    /// to go.
    ///
    /// FIXME: This code is both wrong and ugly, needs rewriting.
    fn find_resonance_phase(&self, f: &netcdf::File) -> Option<f64> {
        use crate::extract::netcdf_fields::*;
        use crate::extract::*;

        let q_values = match extract_1d_array(f, Q) {
            Ok(arr) => arr.to_vec(),
            Err(_) => return None,
        };

        // FIXME: Handle negative and/or non-monotonic q-factors
        let flux_of_q_interp = match &self.flux.state {
            NcFluxState::None => return None, // flux does not exist in the dataset
            _ => {
                match make_interp_type(&self.interp_type)
                    .expect("Has already been build before")
                    .build(&q_values, self.flux.uvalues())
                {
                    Ok(interp) => interp,
                    Err(_) => {
                        eprintln!(
                            "NcHarmonic: q-factor is non-monotonic. Falling back to {:?}",
                            PhaseMethod::fallback()
                        );
                        return None;
                    }
                }
            }
        };

        match flux_of_q_interp.eval(
            &q_values,
            self.flux.uvalues(),
            self.m as f64 / self.n as f64, // q at the resonance
            &mut Accelerator::new(),
        ) {
            Ok(flux_res) => {
                // q=1 returns NaN
                if flux_res.is_finite() {
                    Some(flux_res)
                } else {
                    None
                }
            }
            Err(_) => {
                eprintln!(
                    "NcHarmonic: resonance q is out of bounds. Falling back to {:?}",
                    PhaseMethod::fallback()
                );
                None
            }
        }
    }
}

/// External cache update
impl SingleNcHarmonic {
    /// Updates the interpolated values, since they cannot take place inside the cache without
    /// messing up the whole structure.
    ///
    /// Should always be called together with `cache.update()`, and `update_cache_interps()` should
    /// always be called firsts, since `cache` uses the phase value.
    fn update_cache_interps(&self, flux: f64, cache: &mut NcHarmonicCache) -> Result<()> {
        let acc = cache.flux_acc();

        let new_ampl = match self.alpha_interp.as_ref() {
            Some(i) => i.eval(self.flux.uvalues(), &self.alpha_values, flux, acc)?,
            None => unreachable!("Already checked"),
        };
        let new_dampl = match self.alpha_interp.as_ref() {
            Some(i) => i.eval_deriv(self.flux.uvalues(), &self.alpha_values, flux, acc)?,
            None => unreachable!("Already checked"),
        };
        let new_phase = match self.phase_interp.as_ref() {
            Some(i) => i.eval(self.flux.uvalues(), &self.phase_values, flux, acc)?,
            None => unreachable!("Already checked"),
        };

        cache.set_ampl(new_ampl);
        cache.set_dampl(new_dampl);
        cache.set_phase(new_phase);
        Ok(())
    }

    /// Necessary since cache hits do not call the interpolators, so there is no [`NcFluxState`]
    /// checking.
    fn undefined_evaluation_check(&self, quantity: &str) -> Result<()> {
        if self.flux.state != NcFluxState::Good {
            Err(EqError::UndefinedEvaluation(
                format!("{quantity}({})", self.which).into(),
            ))
        } else {
            Ok(())
        }
    }
}

/// Cache for [`NcHarmonic`] values.
///
/// # Note
///
/// Do not use the same [`NcHarmonicCache`] for both `ψ` and `ψp` evaluations, since this would
/// invalidate the cached values.
///
/// # Example
///
/// ```
/// # use std::path::PathBuf;
/// # use dexter_equilibrium::*;
/// # let path = PathBuf::from("./netcdf.nc");
/// # let harmonic = NcHarmonicBuilder::new(&path, "steffen", 3, 2).build()?;
/// let mut cache = harmonic.get_default_cache();
/// let h = harmonic.h_of_psi(0.0, 0.2, 0.3, 0.0, &mut cache)?;
/// # Ok::<_, EqError>(())
/// ```
///
#[derive(Debug)]
#[non_exhaustive]
pub struct NcHarmonicCache {
    hits: usize,
    misses: usize,
    flux: f64,
    theta: f64,
    zeta: f64,
    ampl: f64,
    dampl: f64, // derivative
    phase: f64,
    modarg: f64,
    cos: f64,
    sin: f64,
    flux_acc: Accelerator,
}

harmonic_cache_counts_getter_impl!(NcHarmonicCache);

impl HarmonicCache for NcHarmonicCache {
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, _: f64) -> bool {
        if (self.flux == flux) && (self.theta == theta) && (self.zeta == zeta) {
            self.hits += 1;
            true
        } else {
            self.misses += 1;
            false
        }
    }

    fn update(&mut self, flux: f64, theta: f64, zeta: f64, _: f64, args: &[f64]) -> Result<()> {
        let m = args.first().expect("Always exists");
        let n = args.get(1).expect("Always exists");

        self.flux = flux;
        self.theta = theta;
        self.zeta = zeta;
        self.modarg = (m * theta - n * zeta + self.phase).rem_euclid(TAU);
        (self.sin, self.cos) = self.modarg.sin_cos();
        Ok(())
    }

    fn ampl(&self) -> f64 {
        self.ampl
    }

    fn sin(&self) -> f64 {
        self.sin
    }

    fn cos(&self) -> f64 {
        self.cos
    }

    fn set_ampl(&mut self, ampl: f64) {
        self.ampl = ampl
    }

    fn set_dampl(&mut self, dampl: f64) {
        self.dampl = dampl
    }

    fn dampl(&self) -> f64 {
        self.dampl
    }

    fn set_phase(&mut self, phase: f64) {
        self.phase = phase
    }

    fn phase(&self) -> f64 {
        self.phase
    }

    fn flux_acc(&mut self) -> &mut Accelerator {
        &mut self.flux_acc
    }
}

impl Default for NcHarmonicCache {
    /// Set all initial values to NaN so first is_updated always returns false.
    fn default() -> Self {
        Self {
            hits: 0,
            misses: 0,
            flux: f64::NAN,
            theta: f64::NAN,
            zeta: f64::NAN,
            modarg: f64::NAN,
            sin: f64::NAN,
            cos: f64::NAN,
            ampl: f64::NAN,
            dampl: f64::NAN,
            phase: f64::NAN,
            flux_acc: Accelerator::new(),
        }
    }
}

/// Intermediate Interpolations. This is effectively where the [`Harmonic`] trait is implemented.
#[rustfmt::skip]
impl SingleNcHarmonic {
    fn ampl(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("α")?;
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t, &self._cache_args)?;
        };
        Ok(cache.ampl())
    }

    fn phase(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("φ")?;
        match self.phase_method {
            PhaseMethod::Zero => Ok(0.0),
            PhaseMethod::Average => Ok(self.phase_average.expect("Exists")),
            PhaseMethod::Resonance => Ok(self.phase_resonance.expect("Exists")),
            PhaseMethod::Custom(phase) => Ok(phase),
            PhaseMethod::Interpolation => {
                if !cache.is_updated(flux, theta, zeta, t) {
                    self.update_cache_interps(flux, cache)?;
                    cache.update(flux, theta, zeta, t, &self._cache_args)?;
                };
                Ok(cache.phase())
            }
        }
    }

    fn h(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("h")?;
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t, &self._cache_args)?;
        };
        Ok(cache.ampl() * cache.cos())
    }

    fn dh_dflux(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("dh")?;
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t, &self._cache_args)?;
        };
        Ok(cache.dampl() * cache.cos())
    }

    fn dh_dtheta(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("dh/dθ")?;
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t, &self._cache_args)?;
        };
        Ok(-self._m * cache.ampl() * cache.sin())
    }

    fn dh_dzeta(&self, flux: f64, theta: f64, zeta: f64, t: f64, cache: &mut NcHarmonicCache) -> Result<f64> {
        self.undefined_evaluation_check("dh/dζ")?;
        if !cache.is_updated(flux, theta, zeta, t) {
            self.update_cache_interps(flux, cache)?;
            cache.update(flux, theta, zeta, t, &self._cache_args)?;
        };
        Ok(self._n * cache.ampl() * cache.sin())
    }
}

impl std::fmt::Debug for SingleNcHarmonic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SingleNcHarmonic")
            .field("NcFlux", &self.flux)
            .field("phase method", &self.phase_method)
            .field("phase average", &self.phase_average)
            .field("phase resonance", &self.phase_resonance)
            .finish()
    }
}

// ===============================================================================================

#[cfg(test)]
impl NcHarmonicCache {
    fn reset(&mut self) {
        *self = Self::default();
    }
}

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
    use is_close::is_close;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn zero_phase_method() {
        use PhaseMethod::Zero;
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            .with_phase_method(Zero)
            .build()
            .unwrap();
        let c = &mut h.get_default_cache();

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
        let c = &mut h.get_default_cache();

        let expected = h.phase_array().mean().unwrap();

        assert!(matches!(h.phase_method(), Average));
        assert!(h.psi_single.phase_average.is_some());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_eq!(h.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), expected);
        assert_eq!(h.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), expected);
    }

    #[test]
    fn resonance_phase_method() {
        use PhaseMethod::Resonance;
        let path = PathBuf::from(POLOIDAL_TEST_NETCDF_PATH);
        let h = NcHarmonicBuilder::new(&path, "steffen", 3, 2)
            .with_phase_method(Resonance)
            .build()
            .unwrap();
        let c = &mut h.get_default_cache();

        assert!(h.psip_single.phase_average.is_none());

        // m=3 and n=2, we expect a resonance at q=n/m=2/3=0.666, which is inbounds in our
        // poloidal_test_netcdf. However, only ψp is a good coordinate.

        // From the create_test_netcdf script:
        let q_res = h.n() as f64 / h.m() as f64;
        let psip_res = TAU * (TAU * q_res).cos(); // we defined ψ = sin(2πψp)
        let expected = h.phase_of_psip(psip_res, 0.1, 0.1, 0.1, c).unwrap();

        assert!(matches!(h.phase_method(), Resonance));
        assert!(h.psip_single.phase_average.is_none());
        assert!(
            h.psip_single
                .phase_resonance
                .is_some_and(|x| is_close!(x, expected))
        );
    }

    #[test]
    fn interpolation_phase_method() {
        use PhaseMethod::Interpolation;
        let h = create_nc_harmonic_builder(POLOIDAL_TEST_NETCDF_PATH)
            .with_phase_method(Interpolation)
            .build()
            .unwrap();
        let c = &mut h.get_default_cache();

        // Calculated with a stable version, on the same dataset
        let expected = 0.8414720460746888;

        assert!(matches!(h.phase_method(), Interpolation));
        assert!(h.psi_single.phase_average.is_none());
        assert!(h.psi_single.phase_resonance.is_none());

        assert!(is_close!(
            h.phase_of_psip(0.1, 0.1, 0.1, 0.0, c).unwrap(),
            expected,
            rel_tol = 1e-5 // unsure why there is this small difference here.
        ));
    }

    #[test]
    fn custom_phase_method() {
        use PhaseMethod::Custom;
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            .with_phase_method(Custom(10.0))
            .build()
            .unwrap();
        let c = &mut h.get_default_cache();

        assert!(matches!(h.phase_method(), Custom(10.0)));
        assert!(h.psi_single.phase_average.is_none());
        assert!(h.psi_single.phase_resonance.is_none());

        assert_eq!(h.phase_of_psi(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
        assert_eq!(h.phase_of_psip(0.01, 0.1, 0.1, 0.0, c).unwrap(), 10.0);
    }

    #[test]
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
        assert_eq!(harmonic.psi_state(), NcFluxState::Good);
        assert_eq!(harmonic.psip_state(), NcFluxState::Bad);
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
        let c = &mut harmonic.get_default_cache();
        assert!(harmonic.ampl_of_psi(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
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
        let c = &mut h.get_default_cache();

        use EqError::UndefinedEvaluation as err;
        assert!(matches!(h.ampl_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.phase_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.h_of_psip(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        c.reset(); // Make sure
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
        assert_eq!(harmonic.psi_state(), NcFluxState::Bad);
        assert_eq!(harmonic.psip_state(), NcFluxState::Good);
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
        let c = &mut harmonic.get_default_cache();
        assert!(harmonic.ampl_of_psip(0.1, 0.1, 0.1, 0.1, c).unwrap().is_finite());
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
        let c = &mut h.get_default_cache();

        use EqError::UndefinedEvaluation as err;
        assert!(matches!(h.ampl_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.phase_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        assert!(matches!(h.h_of_psi(0.1, 0.1, 0.1, 0.1, c), Err(err(..))));
        c.reset(); // Make sure
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
    fn counts() {
        let h = create_nc_harmonic_builder(TEST_NETCDF_PATH)
            // other methods don't touch the cache at all
            .with_phase_method(PhaseMethod::Interpolation)
            .build()
            .unwrap();
        let c = &mut h.get_default_cache();
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
        h.ampl_of_psi(psi, theta, zeta, t, c).unwrap();
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
