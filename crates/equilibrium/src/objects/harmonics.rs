//! Representation of a perturbation 's single harmonic.

use common::array1D_getter_impl;
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::f64::consts::TAU;
use std::path::{Path, PathBuf};

use crate::Harmonic;
use crate::Result;
use crate::cache::HarmonicCache;

/// Defines the calculation method of the phase `φ(ψp)` in an Numerical [`Harmonic`].
#[derive(Default, Debug, Clone)]
pub enum PhaseMethod {
    /// Corresponds to ``φ(ψp) = 0`.
    Zero,
    /// Corresponds to `φ = const = the average of all the values of the phase array`.
    Average,
    /// Corresponds to `φ = const = the value of φ(ψp) at the resonance m/n`.
    #[default]
    Resonance,
    /// Interpolation over the phase array.
    Interpolation,
    /// Use a custom value for φ = const.
    Custom(f64),
}

/// Used to create an [`NcHarmonic`].
#[non_exhaustive]
pub struct NcHarmonicBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
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
    /// # use equilibrium::harmonics;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = harmonics::NcHarmonicBuilder::new(&path, "steffen", 1, 2);
    /// ```
    pub fn new(path: &Path, typ: &str, m: i64, n: i64) -> Self {
        Self {
            path: path.to_path_buf(),
            typ: typ.into(),
            m,
            n,
            phase_method: PhaseMethod::default(),
        }
    }

    /// Creates a new [`NcHarmonic`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::harmonics;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let harmonic = harmonics::NcHarmonicBuilder::new(&path, "cubic", 1, 2).build()?;
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn build(self) -> Result<NcHarmonic> {
        NcHarmonic::build(self)
    }

    /// Sets the phase `φ(ψp)` calculation method.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use equilibrium::harmonics;
    /// # use equilibrium::harmonics::PhaseMethod;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = harmonics::NcHarmonicBuilder::new(&path, "steffen", 1, 2)
    ///     .with_phase_method(PhaseMethod::Interpolation)
    ///     .build()?;
    /// # Ok::<_, equilibrium::EqError>(())
    /// ```
    pub fn with_phase_method(mut self, method: PhaseMethod) -> Self {
        self.phase_method = method;
        self
    }
}

// ===============================================================================================

/// Single perturbation harmonic from a netCDF file.
///
/// The harmonic has the form of `α(ψp) * cos(mθ-nζ+φ(ψp))`, where `α(ψp)` is calculated by
/// interpolation over numerical data, and `φ(ψp)` is calculated as defined by [`PhaseMethod`].
///
/// Should be created with an [`NcHarmonicBuilder`].
#[non_exhaustive]
pub struct NcHarmonic {
    /// Path to the netCDF file.
    path: PathBuf,
    /// 1D [`DynInterpolation`], in case-insensitive string format.
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
    psip_data: Vec<f64>,
    /// The `α` data array.
    alpha_data: Vec<f64>,
    /// The `φ` data array.
    phase_data: Vec<f64>,

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

        let psip_data = extract_1d_array(&f, PSIP_NORM)?.to_vec();
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

        let mut phase_method = builder.phase_method;
        let phase_resonance = match phase_method {
            // Interpolate to find the φ(ψp) value at ψp = m/n.
            PhaseMethod::Resonance => {
                let mut acc = Accelerator::new();
                // Find the ψp where q(ψp) = m/n by creating an inverse q(ψp) -> ψp interpolator
                let resq = (builder.m as f64) / (builder.n as f64);
                let q_data = extract_1d_array::<f64>(&f, Q)?.to_vec();
                let inv_interp = make_interp_type(&builder.typ)?.build(&q_data, &psip_data)?;
                let res = inv_interp.eval(&q_data, &psip_data, resq, &mut acc);
                // If the resonance is outside the wall, use φ=0.
                match res {
                    Ok(psip_res) => Some(
                        phase_interp
                            .eval(&psip_data, &phase_data, psip_res, &mut acc)
                            .unwrap(), // Safe: psip_res is in-bounds
                    ),
                    Err(_) => {
                        phase_method = PhaseMethod::Zero;
                        None
                    }
                }
            }
            _ => None,
        };

        Ok(Self {
            path: path.to_owned(),
            typ: builder.typ,
            _m: builder.m as f64,
            _n: builder.n as f64,
            phase_method,
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
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * cache.cos)
    }

    fn dh_dpsip(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.da_dpsip * cache.cos)
    }

    fn dh_dtheta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * (-self._m) * cache.sin)
    }

    fn dh_dzeta(
        &self,
        psip: f64,
        theta: f64,
        zeta: f64,
        acc: &mut Accelerator,
        cache: &mut HarmonicCache,
    ) -> Result<f64> {
        if !cache.is_updated(psip, theta, zeta) {
            cache.update(self, psip, theta, zeta, acc)?
        };
        Ok(cache.alpha * self._n * cache.sin)
    }

    fn a(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .alpha_interp
            .eval(&self.psip_data, &self.alpha_data, psip, acc)?)
    }

    fn da_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok(self
            .alpha_interp
            .eval_deriv(&self.psip_data, &self.alpha_data, psip, acc)?)
    }

    /// Returns the phase value `φ(ψp)`, depending on the harmonic's [`PhaseMethod`].
    fn phase(&self, psip: f64, acc: &mut Accelerator) -> Result<f64> {
        // Options are always Some when the correct method is set
        match self.phase_method {
            PhaseMethod::Zero => Ok(0.0),
            PhaseMethod::Average => Ok(self.phase_average.expect("is Some")),
            PhaseMethod::Resonance => Ok(self.phase_resonance.expect("is Some")),
            PhaseMethod::Custom(value) => Ok(value),
            PhaseMethod::Interpolation => {
                Ok(self
                    .phase_interp
                    .eval(&self.psip_data, &self.phase_data, psip, acc)?)
            }
        }
    }

    fn mod_arg(&self, psip: f64, theta: f64, zeta: f64, acc: &mut Accelerator) -> Result<f64> {
        Ok((self._m * theta - self._n * zeta + self.phase(psip, acc)?).rem_euclid(TAU))
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
    pub fn phase_average(&self) -> Option<f64> {
        self.phase_average
    }

    /// Returns the [`NcHarmonic`]'s phase value at the resonance φ(ψp) = m/n.
    ///
    /// Returns `None` if the [`NcHarmonic`]'s [`PhaseMethod`] is not [`PhaseMethod::Resonance`].
    pub fn phase_resonance(&self) -> Option<f64> {
        self.phase_resonance
    }

    array1D_getter_impl!(psip_data, psip_data);
    array1D_getter_impl!(a_data, alpha_data);
    array1D_getter_impl!(phase_data, phase_data);
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
