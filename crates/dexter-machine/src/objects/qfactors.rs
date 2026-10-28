//! Representation of a machines's q-factor profile.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_flux, debug_assert_non_negative_psip,
    debug_assert_non_negative_psit, interp_type_getter_impl, netcdf_path_getter_impl,
    netcdf_version_getter_impl,
};
use core::hint::cold_path;
use rsl_interpolation::Accelerator;
use std::path::{Path, PathBuf};

use crate::Interpolation1dType;
use crate::objects::nc_flux::{FluxCoordinateState, NcFlux};
use crate::{EvalError, MachineError};
use crate::{MachineObject, MachineType, Qfactor};
use crate::{MagneticFlux, MagneticFlux::*};
use dexter_common::{DynInterpolator, make_interp};

// ===============================================================================================

/// Analytical q-factor profile of q = 1 and ψt=ψp.
#[non_exhaustive]
pub struct UnityQfactor {
    /// The value of the last closed toroidal flux surface `ψt_last` in Normalized units.
    psit_last: f64,
    /// The value of the last closed poloidal flux surface `ψp_last` in Normalized units.
    psip_last: f64,
}

impl UnityQfactor {
    /// Creates a new `UnityQfactor`.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// // Define ψt=ψt_last=0.45
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let qfactor = UnityQfactor::new(lcfs);
    /// ```
    #[must_use]
    pub fn new(lcfs: MagneticFlux) -> Self {
        Self {
            psit_last: lcfs.value(),
            psip_last: lcfs.value(),
        }
    }
}

impl MachineObject for UnityQfactor {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psit_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }
}

impl Qfactor for UnityQfactor {
    fn psit_last(&self) -> MagneticFlux {
        Toroidal(self.psit_last)
    }

    fn psip_last(&self) -> MagneticFlux {
        Poloidal(self.psip_last)
    }

    fn qlast(&self) -> f64 {
        1.0
    }

    fn qaxis(&self) -> f64 {
        1.0
    }

    fn eval_q(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let last = match flux {
            Toroidal(_) => self.psit_last,
            Poloidal(_) => self.psip_last,
        };
        if flux.value() > last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn eval_q_deriv(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let last = match flux {
            Toroidal(_) => self.psit_last,
            Poloidal(_) => self.psip_last,
        };
        if flux.value() > last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(0.0)
    }

    fn eval_other(
        &self,
        flux: MagneticFlux,
        _: &mut Accelerator,
    ) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let last = match flux {
            Toroidal(_) => self.psit_last,
            Poloidal(_) => self.psip_last,
        };
        let value = flux.value();
        if value > last {
            return Err(EvalError::AnalyticalDomainError);
        };
        match flux {
            Toroidal(_) => Ok(Poloidal(value)),
            Poloidal(_) => Ok(Toroidal(value)),
        }
    }

    fn eval_psit_of_q(&self, _: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        Err(EvalError::UndefinedEvaluation("ψt(q)".into()))
    }

    fn eval_psip_of_q(&self, _: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        Err(EvalError::UndefinedEvaluation("ψp(q)".into()))
    }

    fn eval_deriv_of_other(
        &self,
        flux: MagneticFlux,
        _: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let last = match flux {
            Toroidal(_) => self.psit_last,
            Poloidal(_) => self.psip_last,
        };
        if flux.value() > last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }

    fn eval_deriv_wrt_other(
        &self,
        flux: MagneticFlux,
        _: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let last = match flux {
            Toroidal(_) => self.psit_last,
            Poloidal(_) => self.psip_last,
        };
        if flux.value() > last {
            return Err(EvalError::AnalyticalDomainError);
        };
        Ok(1.0)
    }
}

impl std::fmt::Debug for UnityQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("q-factor of q = 1 and ψt=ψp")
            .field("ψ_last", &self.psit_last)
            .finish()
    }
}

// ===============================================================================================

/// Analytical q-factor of parabolic q(ψt) profile.
///
/// # Note
///
/// The q-factor profile is parabolic *with respect to the toroidal flux `ψt`*.
pub struct ParabolicQfactor {
    /// The `q` value on the magnetic axis.
    qaxis: f64,
    /// The `q` value at the last closed flux surface.
    qlast: f64,
    /// The value of the last closed toroidal flux surface `ψt_last` in Normalized units.
    psit_last: f64,
    /// The value of the last closed poloidal flux surface `ψp_last` in Normalized units.
    psip_last: f64,
}

impl ParabolicQfactor {
    /// Creates a new `ParabolicQfactor`.
    ///
    /// A `ParabolicQfactor` is defined with the help of the [`MagneticFlux`] helper struct,
    /// which changes the position where `qlast` is met.
    ///
    /// # Note
    ///
    /// The q-factor profile is parabolic *with respect to the toroidal flux `ψt`*. The `kind` of `lcfs`
    /// defines the *position* where the value *qlast* is met.
    ///
    /// # Example
    /// ```
    /// # use dexter_machine::*;
    /// // Define q(ψt=ψt_last=0.45) = qlast = 3.8
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    ///
    /// // Define q(ψp=ψp_last=0.19) = qlast = 4.2
    /// let psip_last = MagneticFlux::Poloidal(0.19);
    /// let qfactor = ParabolicQfactor::new(1.1, 4.2, psip_last);
    /// ```
    #[must_use]
    pub fn new(qaxis: f64, qlast: f64, lcfs: MagneticFlux) -> Self {
        // Create two phony qfactors to calculate the other last value
        let psit_last: f64;
        let psip_last: f64;
        match lcfs {
            Toroidal(_psit_last) => {
                let phony_q = Self {
                    qaxis,
                    qlast,
                    psit_last: _psit_last,
                    psip_last: f64::NAN,
                };
                psit_last = phony_q.psit_last;
                psip_last = match phony_q.eval_other(Toroidal(_psit_last), &mut Accelerator::new())
                {
                    Ok(flux) => flux.value(),
                    Err(_) => unreachable!("Analytical formula, cannot fail"),
                };
            }
            Poloidal(_psip_last) => {
                let phony_q = Self {
                    qaxis,
                    qlast,
                    psit_last: f64::NAN,
                    psip_last: _psip_last,
                };
                psip_last = phony_q.psip_last;
                psit_last = phony_q.psit_last_of_psip_last();
            }
        }

        Self {
            qaxis,
            qlast,
            psit_last,
            psip_last,
        }
    }

    /// Helper function to calculate `ψt_last` from `ψp_last` when instantiating `ParabolicQfactor`. This
    /// little maneuver is necessary since `ψt(ψp)` is defined through `ψt_last`, which of course does
    /// not exist yet.
    ///
    /// This formula is derived by solving `q(ψp)` for `ψt_last` and setting `ψp = ψp_last`.
    ///
    /// Should not be used after instantiation.
    #[must_use]
    fn psit_last_of_psip_last(&self) -> f64 {
        let atan_arg = (self.qlast / self.qaxis - 1.0).sqrt();
        let numerator = self.psip_last * (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        let denominator = atan_arg.atan();
        numerator / denominator
    }
}

impl MachineObject for ParabolicQfactor {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psit_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }
}

impl Qfactor for ParabolicQfactor {
    fn psit_last(&self) -> MagneticFlux {
        Toroidal(self.psit_last)
    }

    fn psip_last(&self) -> MagneticFlux {
        Poloidal(self.psip_last)
    }

    fn qaxis(&self) -> f64 {
        self.qaxis
    }

    fn qlast(&self) -> f64 {
        self.qlast
    }

    fn eval_q(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.value()
            > match flux {
                Toroidal(_) => self.psit_last,
                Poloidal(_) => self.psip_last,
            }
        {
            cold_path();
            return Err(EvalError::AnalyticalDomainError);
        };
        let res = match flux {
            Toroidal(psit) => {
                (self.qaxis) + (self.qlast - self.qaxis) * (psit / self.psit_last).powi(2)
            }
            Poloidal(psip) => {
                let tan_arg =
                    (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psit_last;
                self.qaxis + self.qaxis * tan_arg.tan().powi(2)
            }
        };
        Ok(debug_assert_is_finite!(res))
    }

    fn eval_q_deriv(&self, flux: MagneticFlux, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.value()
            > match flux {
                Toroidal(_) => self.psit_last,
                Poloidal(_) => self.psip_last,
            }
        {
            cold_path();
            return Err(EvalError::AnalyticalDomainError);
        };
        let res = match flux {
            Toroidal(psit) => 2.0 * (self.qlast - self.qaxis) * psit / self.psit_last.powi(2),
            Poloidal(psip) => {
                let root = (self.qaxis * (self.qlast - self.qaxis)).sqrt();
                let arg = root * psip / self.psit_last;
                let (sin, cos) = arg.sin_cos();
                let numerator = 2.0 * self.qaxis * root * sin;
                let denominator = self.psit_last * cos.powi(3);
                numerator / denominator
            }
        };
        Ok(debug_assert_is_finite!(res))
    }

    fn eval_other(
        &self,
        flux: MagneticFlux,
        _: &mut Accelerator,
    ) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.value()
            > match flux {
                Toroidal(_) => self.psit_last,
                Poloidal(_) => self.psip_last,
            }
        {
            cold_path();
            return Err(EvalError::AnalyticalDomainError);
        };
        match flux {
            Toroidal(psit) => {
                let atan_arg =
                    psit * (self.qlast - self.qaxis).sqrt() / (self.psit_last * self.qaxis.sqrt());
                let coef = self.psit_last / (self.qaxis * (self.qlast - self.qaxis)).sqrt();
                Ok(Poloidal(debug_assert_is_finite!(coef * atan_arg.atan())))
            }
            Poloidal(psip) => {
                let tan_arg =
                    (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psit_last;
                let coef = self.psit_last * self.qaxis.sqrt() / (self.qlast - self.qaxis).sqrt();
                Ok(Toroidal(debug_assert_is_finite!(coef * tan_arg.tan())))
            }
        }
    }

    fn eval_psit_of_q(&self, q: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        if q > self.qlast {
            return Err(EvalError::AnalyticalDomainError);
        };
        let frac = (q - self.qaxis) / (self.qlast - self.qaxis);
        let psit = self.psit_last * frac.sqrt();
        debug_assert_non_negative_psit!(psit);
        Ok(Toroidal(psit))
    }

    fn eval_psip_of_q(&self, q: f64, _: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        if q > self.qlast {
            return Err(EvalError::AnalyticalDomainError);
        };
        let tan_arg = (q / self.qaxis - 1.0).sqrt();
        let coef = self.psit_last / (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        let psip = coef * tan_arg.atan();
        debug_assert_non_negative_psip!(psip);
        Ok(Poloidal(psip))
    }

    fn eval_deriv_of_other(
        &self,
        flux: MagneticFlux,
        _: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.value()
            > match flux {
                Toroidal(_) => self.psit_last,
                Poloidal(_) => self.psip_last,
            }
        {
            cold_path();
            return Err(EvalError::AnalyticalDomainError);
        };
        match flux {
            Toroidal(psit) => {
                let denom =
                    self.qaxis * self.psit_last.powi(2) + (self.qlast - self.qaxis) * psit.powi(2);
                Ok(debug_assert_is_finite!(self.psit_last.powi(2) / denom))
            }
            Poloidal(psip) => {
                let cos_arg =
                    (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psit_last;
                Ok(debug_assert_is_finite!(self.qaxis / cos_arg.cos().powi(2)))
            }
        }
    }

    fn eval_deriv_wrt_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let other = self.eval_other(flux, acc)?;
        self.eval_deriv_of_other(other, &mut Accelerator::new())
    }
}

impl std::fmt::Debug for ParabolicQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParabolicQfactor: q-factor of parabolic q(ψt) profile.")
            .field("ψt_last", &self.psit_last)
            .field("ψp_last", &self.psip_last)
            .field("qaxis", &self.qaxis)
            .field("qlast", &self.qlast)
            .finish()
    }
}

// ===============================================================================================

/// Used to create a [`NcQfactor`].
///
/// Exists for future configuration flexibility.
#[derive(Debug)]
pub struct NcQfactorBuilder {
    /// Path to the netCDF file.
    path: PathBuf,
    /// The interpolation type.
    interp_type: Interpolation1dType,
}

impl NcQfactorBuilder {
    /// Creates a new `NcQfactorBuilder` from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, Interpolation1dType::Cubic);
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: Interpolation1dType) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type,
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_machine::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, Interpolation1dType::Akima).build()?;
    /// # Ok::<_, MachineError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`MachineError`] if it fails to build the [`NcQfactor`].
    pub fn build(self) -> Result<NcQfactor, MachineError> {
        NcQfactor::build(self)
    }
}

/// Numerical q-factor profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcQfactorBuilder`].
///
/// If either `psit_norm` or `psip_norm` is missing from the netCDF file, it is calculated from the
/// other by integrating `q(ψp)` or `ι(ψt)` respectively. In the case that the calculated values are
/// monotonic, the other flux can be used as a flux coordinate as well.
pub struct NcQfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The interpolation type.
    interp_type: Interpolation1dType,

    /// The toroidal flux coordinate.
    psit: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// `ψp(ψt)` interpolator.
    psip_of_psit_interp: Option<DynInterpolator>,
    /// `ψt(ψp)` interpolator.
    psit_of_psip_interp: Option<DynInterpolator>,

    /// The `q` values.
    q_values: Box<[f64]>,
    /// `q(ψt)` interpolator.
    q_of_psit_interp: Option<DynInterpolator>,
    /// `q(ψp)` interpolator.
    q_of_psip_interp: Option<DynInterpolator>,

    /// `ψt(q)` interpolator.
    psit_of_q_interp: Option<DynInterpolator>,
    /// `ψp(q)` interpolator.
    psip_of_q_interp: Option<DynInterpolator>,
}

/// Creation.
impl NcQfactor {
    /// Constructs an `NcQfactor` from an [`NcQfactorBuilder`].
    fn build(builder: NcQfactorBuilder) -> Result<Self, MachineError> {
        use crate::extract;
        use crate::extract::netcdf_fields::NC_Q;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let mut psit = NcFlux::toroidal(&file);
        let mut psip = NcFlux::poloidal(&file);
        let q_values = extract::array_1d(&file, NC_Q)?.to_vec();

        // If no `ψt` values found, integrate `q(ψp)` and create them. We create a temporary
        // `q_of_psip_interp` since we cannot correctly define it yet.
        //
        // This operation cannot fail. If one of the fluxes is missing, then the other one is
        // guaranteed to be `Good` and therefore the interpolation and integration can be defined.
        if psit.state() == FluxCoordinateState::NoValues {
            let acc = &mut Accelerator::new();
            let psip_values = psip.values().expect("At least one of the fluxes exists");
            let q_of_psip_interp = make_interp(builder.interp_type, psip_values, &q_values)?;
            let psit_values: Vec<f64> = psip_values
                .iter()
                .map(|psip_value| {
                    q_of_psip_interp
                        .eval_integ(psip_values, &q_values, 0.0, *psip_value, acc)
                        .expect("Cannot fail by definition of the interpolator")
                })
                .collect();
            psit = NcFlux::from_raw_values(&psit_values);
        }

        // If no `ψp` values found, integrate `ι(ψt)` and create them. For this we must create a
        // temporary `i_of_psit` interpolator.
        //
        // This operation cannot fail. If one of the fluxes is missing, then the other one is
        // guaranteed to be `Good` and therefore the interpolation and integration can be defined.
        if psip.state() == FluxCoordinateState::NoValues {
            let acc = &mut Accelerator::new();
            let psit_values = psit.values().expect("At least one of the fluxes exists");
            let i_values: Vec<f64> = q_values.iter().map(|q| q.recip()).collect();
            let i_of_psit_interp = make_interp(builder.interp_type, psit_values, &i_values)?;
            let psip_values: Vec<f64> = psit_values
                .iter()
                .map(|psit_value| {
                    i_of_psit_interp
                        .eval_integ(psit_values, &i_values, 0.0, *psit_value, acc)
                        .expect("Cannot fail by definition of the interpolator")
                })
                .collect();
            psip = NcFlux::from_raw_values(&psip_values);
        }

        // Create interpolators, if possible
        use FluxCoordinateState::Good;
        let psip_of_psit_interp =
            if (psit.state() == Good) & (psip.state() != FluxCoordinateState::NoValues) {
                Some(make_interp(
                    builder.interp_type,
                    psit.uvalues(),
                    psip.uvalues(),
                )?)
            } else {
                None
            };
        let psit_of_psip_interp =
            if (psip.state() == Good) & (psit.state() != FluxCoordinateState::NoValues) {
                Some(make_interp(
                    builder.interp_type,
                    psip.uvalues(),
                    psit.uvalues(),
                )?)
            } else {
                None
            };

        let q_of_psit_interp = match psit.state() {
            Good => Some(make_interp(builder.interp_type, psit.uvalues(), &q_values)?),
            _ => None,
        };
        let q_of_psip_interp = match psip.state() {
            Good => Some(make_interp(builder.interp_type, psip.uvalues(), &q_values)?),
            _ => None,
        };

        // If flux values exist, we must also check if q is monotonic
        let psit_of_q_interp = match psit.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp_type, &q_values, psit.uvalues()).ok(),
        };
        let psip_of_q_interp = match psip.state() {
            FluxCoordinateState::NoValues => None,
            _ => make_interp(builder.interp_type, &q_values, psip.uvalues()).ok(),
        };

        Ok(Self {
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            psit,
            psip,
            psip_of_psit_interp,
            psit_of_psip_interp,
            q_values: q_values.into_boxed_slice(),
            q_of_psit_interp,
            q_of_psip_interp,
            psit_of_q_interp,
            psip_of_q_interp,
        })
    }
}

impl MachineObject for NcQfactor {
    fn machine_type(&self) -> MachineType {
        MachineType::Numerical
    }

    fn psit_state(&self) -> FluxCoordinateState {
        self.psit.state()
    }

    fn psip_state(&self) -> FluxCoordinateState {
        self.psip.state()
    }
}

impl Qfactor for NcQfactor {
    fn psit_last(&self) -> MagneticFlux {
        Toroidal(
            self.psit
                .last_value()
                .expect("'psit' values have been calculated at this point"),
        )
    }

    fn psip_last(&self) -> MagneticFlux {
        Poloidal(
            self.psip
                .last_value()
                .expect("'psip' values have been calculated at this point"),
        )
    }

    fn qlast(&self) -> f64 {
        self.q_values
            .last()
            .copied()
            .expect("'q' array is non-empty")
    }

    fn qaxis(&self) -> f64 {
        self.q_values
            .first()
            .copied()
            .expect("'q' array is non-empty")
    }

    fn eval_q(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (val, xa, interp) = match flux {
            Toroidal(val) => (val, self.psit.uvalues(), self.q_of_psit_interp.as_ref()),
            Poloidal(val) => (val, self.psip.uvalues(), self.q_of_psip_interp.as_ref()),
        };
        let ya = &self.q_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(interp.eval(xa, ya, val, acc)?))
        } else {
            cold_path();
            let msg = format!("q({})", flux.symbol());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }

    fn eval_q_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let (val, xa, interp) = match flux {
            Toroidal(val) => (val, self.psit.uvalues(), self.q_of_psit_interp.as_ref()),
            Poloidal(val) => (val, self.psip.uvalues(), self.q_of_psip_interp.as_ref()),
        };
        let ya = &self.q_values;
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(
                interp.eval_deriv(xa, ya, val, acc)?
            ))
        } else {
            cold_path();
            let msg = format!("dq({})/d({})", flux.symbol(), flux.symbol());
            Err(EvalError::UndefinedEvaluation(msg))
        }
    }

    fn eval_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<MagneticFlux, EvalError> {
        debug_assert_non_negative_flux!(flux);
        match flux {
            Toroidal(_) if self.psip_of_psit_interp.is_none() => {
                cold_path();
                return Err(EvalError::UndefinedEvaluation("ψp(ψt)".into()));
            }
            Poloidal(_) if self.psit_of_psip_interp.is_none() => {
                cold_path();
                return Err(EvalError::UndefinedEvaluation("ψt(ψp)".into()));
            }
            _ => (),
        }
        // This cannot panic. The above check ensures the flux values exist.
        #[rustfmt::skip]
        let (val, xa, ya, interp) = match flux {
            Toroidal(val) => (val, self.psit.uvalues(), self.psip.uvalues(), self.psip_of_psit_interp.as_ref()),
            Poloidal(val) => (val, self.psip.uvalues(), self.psit.uvalues(), self.psit_of_psip_interp.as_ref()),
        };
        let value = if let Some(interp) = interp {
            debug_assert_is_finite!(interp.eval(xa, ya, val, acc)?)
        } else {
            let msg = match flux {
                Toroidal(_) => "ψp(ψt)",
                Poloidal(_) => "ψt(ψp)",
            };
            cold_path();
            return Err(EvalError::UndefinedEvaluation(msg.into()));
        };
        match flux {
            Toroidal(_) => Ok(Poloidal(value)),
            Poloidal(_) => Ok(Toroidal(value)),
        }
    }

    fn eval_psit_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        match self.psit_of_q_interp.as_ref() {
            Some(interp) => Ok(Toroidal(debug_assert_is_finite!(interp.eval(
                &self.q_values,
                self.psit.uvalues(),
                q,
                acc
            )?))),
            None => Err(EvalError::UndefinedEvaluation("ψt(q)".into())),
        }
    }

    fn eval_psip_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError> {
        match self.psip_of_q_interp.as_ref() {
            Some(interp) => Ok(Poloidal(debug_assert_is_finite!(interp.eval(
                &self.q_values,
                self.psip.uvalues(),
                q,
                acc
            )?))),
            None => Err(EvalError::UndefinedEvaluation("ψp(q)".into())),
        }
    }

    fn eval_deriv_of_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        match flux {
            Toroidal(_) if self.psit_of_psip_interp.is_none() => {
                return Err(EvalError::UndefinedEvaluation("ψp(ψt)".into()));
            }
            Poloidal(_) if self.psip_of_psit_interp.is_none() => {
                return Err(EvalError::UndefinedEvaluation("ψt(ψp)".into()));
            }
            _ => (),
        }
        // This cannot panic. The above check ensures the flux values exist.
        #[rustfmt::skip]
        let (val, xa, ya, interp) = match flux {
            Toroidal(val) => (val, self.psit.uvalues(), self.psip.uvalues(), self.psip_of_psit_interp.as_ref()),
            Poloidal(val) => (val, self.psip.uvalues(), self.psit.uvalues(), self.psit_of_psip_interp.as_ref()),
        };
        if let Some(interp) = interp {
            Ok(debug_assert_is_finite!(
                interp.eval_deriv(xa, ya, val, acc)?
            ))
        } else {
            let msg = match flux {
                Toroidal(_) => "ψp(ψt)",
                Poloidal(_) => "ψt(ψp)",
            };
            Err(EvalError::UndefinedEvaluation(msg.into()))
        }
    }

    fn eval_deriv_wrt_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        let other = self.eval_other(flux, acc)?;
        self.eval_deriv_of_other(other, &mut Accelerator::new())
    }
}

/// Getters.
impl NcQfactor {
    /// Returns a reference to the toroidal flux's values.
    #[must_use]
    #[expect(clippy::missing_panics_doc, reason = "expect")]
    pub fn psit_values(&self) -> &[f64] {
        self.psit.values().expect("always exist by construction")
    }

    /// Returns a reference to the poloidal flux's values.
    #[must_use]
    #[expect(clippy::missing_panics_doc, reason = "expect")]
    pub fn psip_values(&self) -> &[f64] {
        self.psip.values().expect("always exist by construction")
    }

    /// Returns a reference to the `q` values.
    #[must_use]
    pub fn q_values(&self) -> &[f64] {
        &self.q_values
    }

    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    interp_type_getter_impl!(1);
}

impl std::fmt::Debug for NcQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcQfactor")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("interpolation type", &self.interp_type())
            .field("ψt", &self.psit)
            .field("ψp", &self.psip)
            .field("qaxis", &self.qaxis())
            .field("qlast", &self.qlast())
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::*;

    pub(super) fn create_nc_qfactor(path_str: &str) -> NcQfactor {
        let path = PathBuf::from(&path_str);
        let builder = NcQfactorBuilder::new(&path, Interpolation1dType::Akima);
        builder.build().unwrap()
    }

    /// Make sure that dψt(ψp)/dψp and q(ψt) are close enough.
    pub(super) fn test_dpsit_dpsip_q_closeness<Q: Qfactor>(qfactor: &Q) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psips = Array1::linspace(
            0.02 * qfactor.psip_last().value(),
            0.98 * qfactor.psip_last().value(),
            100,
        );

        let mut acc = Accelerator::new();
        for psip in psips.iter().copied() {
            let psip = Poloidal(psip);
            let q = qfactor.eval_q(psip, &mut acc).unwrap();
            let dpsit_dpsip = qfactor.eval_deriv_of_other(psip, &mut acc).unwrap();
            assert_relative_eq!(q, dpsit_dpsip, epsilon = qfactor.qlast() * 1e-4)
        }
    }

    /// Make sure that dψp(ψt)/dψt and i(ψt) are close enough.
    pub(super) fn test_dpsip_dpsit_iota_closeness<Q: Qfactor>(qfactor: &Q) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psits = Array1::linspace(
            0.02 * qfactor.psit_last().value(),
            0.98 * qfactor.psit_last().value(),
            100,
        );

        let mut acc = Accelerator::new();
        for psit in psits.iter().copied() {
            let psit = Toroidal(psit);
            let iota = qfactor.eval_iota(psit, &mut acc).unwrap();
            let dpsip_dpsit = qfactor.eval_deriv_of_other(psit, &mut acc).unwrap();
            assert_relative_eq!(iota, dpsip_dpsit, epsilon = qfactor.qlast() * 1e-4)
        }
    }

    pub(super) fn test_eval_deriv_methods<Q: Qfactor>(qfactor: &Q) {
        assert!(qfactor.psit_state() == FluxCoordinateState::Good);
        assert!(qfactor.psip_state() == FluxCoordinateState::Good);

        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psits = Array1::linspace(
            0.02 * qfactor.psit_last().value(),
            0.98 * qfactor.psit_last().value(),
            100,
        );

        let mut acc = Accelerator::new();
        for psit in psits.iter().copied() {
            let psit = Toroidal(psit);
            let psip = qfactor.eval_other(psit, &mut acc).unwrap();
            let dpsip_dpsit_1 = qfactor.eval_deriv_of_other(psit, &mut acc).unwrap();
            let dpsip_dpsit_2 = qfactor.eval_deriv_wrt_other(psip, &mut acc).unwrap();
            assert_relative_eq!(
                dpsip_dpsit_1,
                dpsip_dpsit_2,
                epsilon = qfactor.qlast() * 1e-4
            )
        }

        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psips = Array1::linspace(
            0.02 * qfactor.psip_last().value(),
            0.98 * qfactor.psip_last().value(),
            100,
        );

        let mut acc = Accelerator::new();
        for psip in psips.iter().copied() {
            let psip = Toroidal(psip);
            let psit = qfactor.eval_other(psip, &mut acc).unwrap();
            let dpsit_dpsip_1 = qfactor.eval_deriv_of_other(psip, &mut acc).unwrap();
            let dpsit_dpsip_2 = qfactor.eval_deriv_wrt_other(psit, &mut acc).unwrap();
            assert_relative_eq!(
                dpsit_dpsip_1,
                dpsit_dpsip_2,
                epsilon = qfactor.qlast() * 1e-4
            )
        }
    }
}

#[cfg(test)]
mod test_unity_qfactor {
    use crate::objects::qfactors::test_utils::test_eval_deriv_methods;

    use super::*;

    #[test]
    fn domain_error() {
        let lcfs = Toroidal(0.1);
        let qfactor = UnityQfactor::new(lcfs);

        let acc = &mut Accelerator::new();
        assert!(qfactor.eval_q(Toroidal(0.2), acc).is_err());
        assert!(qfactor.eval_q(Poloidal(0.2), acc).is_err());
        test_eval_deriv_methods(&qfactor);
    }
}

#[cfg(test)]
mod test_parabolic_qfactor {
    use crate::{MagneticFluxKind, objects::qfactors::test_utils::test_eval_deriv_methods};

    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psit_last() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Toroidal(0.45));
        assert_abs_diff_eq!(qfactor.psit_last().value(), 0.45);
        assert_relative_eq!(
            qfactor.psip_last().value(),
            0.25921022097041035,
            epsilon = 1e-12
        );

        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
    }

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psip_last() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Poloidal(0.25921022097041035));
        assert_abs_diff_eq!(qfactor.psip_last().value(), 0.25921022097041035);
        assert_relative_eq!(qfactor.psit_last().value(), 0.45, epsilon = 1e-12);

        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
    }

    #[test]
    fn inverse_q() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Toroidal(0.45));
        let acc = &mut Accelerator::new();

        let flux = Toroidal(0.01);
        let q_of_psit = qfactor.eval_q(flux, acc).unwrap();
        let psit_inverse = qfactor.eval_psit_of_q(q_of_psit, acc).unwrap();
        assert_eq!(flux.kind(), MagneticFluxKind::Toroidal);
        assert_relative_eq!(flux.value(), psit_inverse.value(), epsilon = 1e-12);

        let flux = Poloidal(0.1);
        let q_of_psip = qfactor.eval_q(flux, acc).unwrap();
        let psip_inverse = qfactor.eval_psip_of_q(q_of_psip, acc).unwrap();
        assert_eq!(flux.kind(), MagneticFluxKind::Poloidal);
        assert_relative_eq!(flux.value(), psip_inverse.value(), epsilon = 1e-12);
    }

    #[test]
    fn eval_derivative_methods() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Toroidal(0.45));
        test_eval_deriv_methods(&qfactor);
    }
}

#[cfg(test)]
mod test_derivatives_closeness {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn unity_qfactor_dpsit_dpsip_q_closeness() {
        let qfactor = UnityQfactor::new(Poloidal(0.45));
        test_dpsit_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn unity_qfactor_dpsip_dpsit_iota_closeness() {
        let qfactor = UnityQfactor::new(Toroidal(0.45));
        test_dpsip_dpsit_iota_closeness(&qfactor);
    }

    #[test]
    fn parabolic_qfactor_dpsit_dpsip_q_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Toroidal(0.45));
        test_dpsit_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn parabolic_qfactor_dpsip_dpsit_iota_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, Toroidal(0.45));
        test_dpsip_dpsit_iota_closeness(&qfactor);
    }

    #[test]
    fn nc_qfactor_dpsit_dpsip_q_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        test_dpsit_dpsip_q_closeness(&qfactor);
    }

    #[test]
    fn nc_qfactor_dpsip_dpsit_iota_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        test_dpsip_dpsit_iota_closeness(&qfactor);
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::{TEST_NETCDF_PATH, TOROIDAL_TEST_NETCDF_PATH};

    use super::test_utils::*;
    use super::*;

    #[test]
    fn eval_derivative_methods() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        test_eval_deriv_methods(&qfactor);
    }

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Bad);

        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Good);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Bad);
        assert!(qfactor.q_of_psit_interp.is_some());
        assert!(qfactor.q_of_psip_interp.is_none());
        assert!(qfactor.psip_of_psit_interp.is_some());
        assert!(qfactor.psit_of_psip_interp.is_none());

        let _ = qfactor.psit_values();
        let _ = qfactor.psip_values();
    }

    #[test]
    fn good_psit_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Toroidal(0.01);
        assert!(qfactor.eval_q(flux, &mut acc).unwrap().is_finite());
        assert!(qfactor.eval_iota(flux, &mut acc).unwrap().is_finite());
        assert!(
            qfactor
                .eval_other(flux, &mut acc)
                .unwrap()
                .value()
                .is_finite()
        );
        // Need "good" poloidal flux
        assert!(
            qfactor
                .eval_deriv_of_other(Toroidal(0.01), &mut acc)
                .is_err()
        );
        assert!(
            qfactor
                .eval_deriv_of_other(Poloidal(0.01), &mut acc)
                .is_err()
        );
    }

    #[test]
    fn bad_psip_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Poloidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.eval_q(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_iota(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_other(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_deriv_of_other(flux, &mut acc), Err(err(..)));
    }
}

#[cfg(test)]
mod test_poloidal_nc_evals {
    use crate::extract::POLOIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Bad);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);

        assert_eq!(qfactor.psit_state(), FluxCoordinateState::Bad);
        assert_eq!(qfactor.psip_state(), FluxCoordinateState::Good);
        assert!(qfactor.q_of_psit_interp.is_none());
        assert!(qfactor.q_of_psip_interp.is_some());
        assert!(qfactor.psip_of_psit_interp.is_none());
        assert!(qfactor.psit_of_psip_interp.is_some());

        let _ = qfactor.psit_values();
        let _ = qfactor.psip_values();
    }

    #[test]
    fn good_psip_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Poloidal(0.01);
        assert!(qfactor.eval_q(flux, &mut acc).unwrap().is_finite());
        assert!(qfactor.eval_iota(flux, &mut acc).unwrap().is_finite());
        assert!(
            qfactor
                .eval_other(flux, &mut acc)
                .unwrap()
                .value()
                .is_finite()
        );
        // Need "good" toroidal flux
        assert!(
            qfactor
                .eval_deriv_of_other(Toroidal(0.01), &mut acc)
                .is_err()
        );
        assert!(
            qfactor
                .eval_deriv_of_other(Poloidal(0.01), &mut acc)
                .is_err()
        );
    }

    #[test]
    fn bad_psit_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        let flux = Toroidal(0.01);
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.eval_q(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_iota(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_other(flux, &mut acc), Err(err(..)));
        matches!(qfactor.eval_deriv_of_other(flux, &mut acc), Err(err(..)));
    }
}
