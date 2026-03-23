//! Representation of an equilibrium's q-factor profile.

use crate::{
    debug_assert_is_finite, debug_assert_non_negative_psi, debug_assert_non_negative_psip,
    equilibrium_type_getter_impl, fluxes_state_getter_impl, fluxes_values_array_getter_impl,
    interp_type_getter_impl, lcfs_getter_impl, netcdf_path_getter_impl, netcdf_version_getter_impl,
};
use dexter_common::vec_to_array1D_getter_impl;
use ndarray::Array1;
use rsl_interpolation::{Accelerator, DynInterpolation, InterpType, make_interp_type};
use std::path::{Path, PathBuf};

use super::debug_assert_all_finite_values;
use crate::objects::nc_flux::{NcFlux, NcFluxState};
use crate::{EqError, EvalError};
use crate::{EquilibriumType, FluxCommute, LastClosedFluxSurface, Qfactor};

// ===============================================================================================

/// Analytical q-factor profile of q = 1 and ψ=ψp.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[non_exhaustive]
pub struct UnityQfactor {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
}

impl UnityQfactor {
    /// Creates a new `UnityQfactor`.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// let qfactor = UnityQfactor::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self {
            equilibrium_type: EquilibriumType::Analytical,
        }
    }

    equilibrium_type_getter_impl!();
}

impl FluxCommute for UnityQfactor {
    fn psip_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(psi)
    }

    fn psi_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(psip)
    }
}

impl Qfactor for UnityQfactor {
    fn q_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(1.0)
    }

    fn q_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(1.0)
    }

    fn dpsip_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(1.0)
    }

    fn dpsi_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        Ok(1.0)
    }
}

impl std::fmt::Debug for UnityQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("q-factor of q = 1 and ψ=ψp").finish()
    }
}

// ===============================================================================================

/// Analytical q-factor of parabolic q(ψ) profile.
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
pub struct ParabolicQfactor {
    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The `q` value on the magnetic axis.
    qaxis: f64,
    /// The `q` value at the last closed flux surface.
    qlast: f64,
    /// The value of the last closed toroidal flux surface `ψ_last` in Normalized units.
    psi_last: f64,
    /// The value of the last closed poloidal flux surface `ψp_last` in Normalized units.
    psip_last: f64,
}

impl ParabolicQfactor {
    /// Creates a new `ParabolicQfactor`.
    ///
    /// A `ParabolicQfactor` is defined with the help of the [`LastClosedFluxSurface`] helper struct,
    /// which changes the position where `qlast` is met.
    ///
    /// # Example
    /// ```
    /// # use dexter_equilibrium::*;
    /// // Define q(ψ=ψ_last=0.45) = qlast = 3.8
    /// let psi_last = LastClosedFluxSurface::Toroidal(0.45);
    /// let qfactor = ParabolicQfactor::new(1.1, 3.8, psi_last);
    ///
    /// // Define q(ψp=ψp_last=0.19) = qlast = 4.2
    /// let psip_last = LastClosedFluxSurface::Poloidal(0.19);
    /// let qfactor = ParabolicQfactor::new(1.1, 4.2, psip_last);
    /// ```
    #[must_use]
    pub fn new(qaxis: f64, qlast: f64, lcfs: LastClosedFluxSurface) -> Self {
        // Create two phony qfactors to calculate the other last value
        let psi_last: f64;
        let psip_last: f64;
        match lcfs {
            LastClosedFluxSurface::Toroidal(_psi_last) => {
                let phony_q = Self {
                    equilibrium_type: EquilibriumType::Analytical,
                    qaxis,
                    qlast,
                    psi_last: _psi_last,
                    psip_last: f64::NAN,
                };
                psi_last = phony_q.psi_last;
                psip_last = match phony_q.psip_of_psi(_psi_last, &mut Accelerator::new()) {
                    Ok(value) => value,
                    Err(_) => unreachable!("Analytical formula, cannot fail"),
                };
            }
            LastClosedFluxSurface::Poloidal(_psip_last) => {
                let phony_q = Self {
                    equilibrium_type: EquilibriumType::Analytical,
                    qaxis,
                    qlast,
                    psi_last: f64::NAN,
                    psip_last: _psip_last,
                };
                psip_last = phony_q.psip_last;
                psi_last = phony_q.psi_last_of_psip_last();
            }
        }

        Self {
            equilibrium_type: EquilibriumType::Analytical,
            qaxis,
            qlast,
            psi_last,
            psip_last,
        }
    }

    equilibrium_type_getter_impl!();

    /// Returns the value of `q` on the magnetic axis.
    #[must_use]
    pub fn qaxis(&self) -> f64 {
        self.qaxis
    }

    /// Returns the value of `q` at the last closed flux surface.
    #[must_use]
    pub fn qlast(&self) -> f64 {
        self.qlast
    }

    /// Returns the value of the last closed toroidal flux surface `ψ_last`.
    #[must_use]
    pub fn psi_last(&self) -> f64 {
        self.psi_last
    }

    /// Returns the value of the last closed poloidal flux surface `ψp_last`.
    #[must_use]
    pub fn psip_last(&self) -> f64 {
        self.psip_last
    }

    /// Helper function to calculate `ψ_last` from `ψp_last` when instantiating `ParabolicQfactor`. This
    /// little maneuver is necessary since `ψ(ψp)` is defined through `ψ_last`, which of course does
    /// not exist yet.
    ///
    /// This formula is derived by solving `q(ψp)` for `ψ_last` and setting `ψp = ψp_last`.
    ///
    /// Should not be used after instantiation.
    #[must_use]
    fn psi_last_of_psip_last(&self) -> f64 {
        let atan_arg = (self.qlast / self.qaxis - 1.0).sqrt();
        let numerator = self.psip_last * (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        let denominator = atan_arg.atan();
        numerator / denominator
    }
}

impl FluxCommute for ParabolicQfactor {
    fn psip_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        let atan_arg = psi * (self.qlast - self.qaxis).sqrt() / (self.psi_last * self.qaxis.sqrt());
        let coef = self.psi_last / (self.qaxis * (self.qlast - self.qaxis)).sqrt();
        Ok(debug_assert_is_finite!(coef * atan_arg.atan()))
    }

    fn psi_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        let tan_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        let coef = self.psi_last * self.qaxis.sqrt() / (self.qlast - self.qaxis).sqrt();
        Ok(debug_assert_is_finite!(coef * tan_arg.tan()))
    }
}

// TODO: Cache reoccurring values when sure the formulas are correct.
impl Qfactor for ParabolicQfactor {
    fn q_of_psi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        Ok(debug_assert_is_finite!(
            (self.qaxis) + (self.qlast - self.qaxis) * (psi / self.psi_last).powi(2)
        ))
    }

    fn q_of_psip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        let tan_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        Ok(debug_assert_is_finite!(
            self.qaxis + self.qaxis * tan_arg.tan().powi(2)
        ))
    }

    fn dpsip_dpsi(&self, psi: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        let denom = self.qaxis * self.psi_last.powi(2) + (self.qlast - self.qaxis) * psi.powi(2);
        Ok(debug_assert_is_finite!(self.psi_last.powi(2) / denom))
    }

    fn dpsi_dpsip(&self, psip: f64, _: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        let cos_arg = (self.qaxis * (self.qlast - self.qaxis)).sqrt() * psip / self.psi_last;
        Ok(debug_assert_is_finite!(self.qaxis / cos_arg.cos().powi(2)))
    }
}

impl std::fmt::Debug for ParabolicQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParabolicQfactor: q-factor of parabolic q(ψ) profile.")
            .field("equilibrium_type", &self.equilibrium_type)
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
    /// 1D [`DynInterpolation`] type (case-insensitive).
    interp_type: String,
}

impl NcQfactorBuilder {
    /// Creates a new [`NcQfactorBuilder`] from a netCDF file at `path`, with `interp_type`
    /// interpolation type.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let builder = NcQfactorBuilder::new(&path, "cubic");
    /// ```
    #[must_use]
    pub fn new(path: &Path, interp_type: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            interp_type: interp_type.into(),
        }
    }

    /// Creates a new [`NcQfactor`] with the Builder's configuration.
    ///
    /// # Example
    /// ```
    /// # use std::path::PathBuf;
    /// # use dexter_equilibrium::*;
    /// let path = PathBuf::from("./netcdf.nc");
    /// let qfactor = NcQfactorBuilder::new(&path, "cubic").build()?;
    /// # Ok::<_, EqError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`EqError`] if it fails to build the [`NcQfactor`].
    pub fn build(self) -> Result<NcQfactor, EqError> {
        NcQfactor::build(self)
    }
}

/// Numerical q-factor profile reconstructed from a netCDF file.
///
/// Related quantities are computed by interpolating over the data arrays.
///
/// Should be created with an [`NcQfactorBuilder`].
pub struct NcQfactor {
    /// Path to the netCDF file.
    path: PathBuf,
    /// netCDF's [`semver::Version`].
    netcdf_version: semver::Version,

    /// The object's equilibrium type.
    equilibrium_type: EquilibriumType,
    /// The interpolation type.
    interp_type: String,

    /// The toroidal flux coordinate.
    psi: NcFlux,
    /// The poloidal flux coordinate.
    psip: NcFlux,

    /// `ψp(ψ)` interpolator.
    psip_of_psi_interp: Option<DynInterpolation<f64>>,
    /// `ψ(ψp)` interpolator.
    psi_of_psip_interp: Option<DynInterpolation<f64>>,

    /// The `q` values.
    q_values: Vec<f64>,
    /// `q(ψ)` interpolator.
    q_of_psi_interp: Option<DynInterpolation<f64>>,
    /// `q(ψp)` interpolator.
    q_of_psip_interp: Option<DynInterpolation<f64>>,
}

/// Creation.
impl NcQfactor {
    /// Constructs an [`NcQfactor`] from an [`NcQfactorBuilder`].
    fn build(builder: NcQfactorBuilder) -> Result<Self, EqError> {
        use crate::extract;
        use crate::extract::netcdf_fields::NC_Q;

        // Make path absolute for display purposes.
        let path = std::path::absolute(builder.path)?;
        let file = extract::open(&path)?;
        let netcdf_version = extract::version(&file)?;

        let psi = NcFlux::toroidal(&file);
        let psip = NcFlux::poloidal(&file);
        let q_values = extract::array_1d(&file, NC_Q)?.to_vec();

        debug_assert_all_finite_values(&q_values);

        // TODO: Integrate q/ι in case one of the fluxes is missing.

        // Create interpolators, if possible
        use NcFluxState::Good;
        let psip_of_psi_interp = if (psi.state() == Good) & (psip.state() != NcFluxState::None) {
            Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), psip.uvalues())?)
        } else {
            None
        };
        let psi_of_psip_interp = if (psip.state() == Good) & (psi.state() != NcFluxState::None) {
            Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), psi.uvalues())?)
        } else {
            None
        };

        let q_of_psi_interp = match psi.state() {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psi.uvalues(), &q_values)?),
            _ => None,
        };
        let q_of_psip_interp = match psip.state() {
            Good => Some(make_interp_type(&builder.interp_type)?.build(psip.uvalues(), &q_values)?),
            _ => None,
        };

        Ok(Self {
            equilibrium_type: EquilibriumType::Numerical,
            netcdf_version,
            path,
            interp_type: builder.interp_type,
            psi,
            psip,
            psip_of_psi_interp,
            psi_of_psip_interp,
            q_values,
            q_of_psi_interp,
            q_of_psip_interp,
        })
    }
}

impl FluxCommute for NcQfactor {
    fn psip_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.psip_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                self.psip.uvalues(),
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψp(ψ)".into())),
        }
    }

    fn psi_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.psi_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                self.psi.uvalues(),
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("ψ(ψp)".into())),
        }
    }
}

impl Qfactor for NcQfactor {
    fn q_of_psi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.q_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psi.uvalues(),
                &self.q_values,
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("q(ψ)".into())),
        }
    }

    fn q_of_psip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.q_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval(
                self.psip.uvalues(),
                &self.q_values,
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("q(ψp)".into())),
        }
    }

    fn dpsip_dpsi(&self, psi: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psi!(psi);
        match self.psip_of_psi_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psi.uvalues(),
                self.psip.uvalues(),
                psi,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dψp(ψ)/dψ".into())),
        }
    }

    fn dpsi_dpsip(&self, psip: f64, acc: &mut Accelerator) -> Result<f64, EvalError> {
        debug_assert_non_negative_psip!(psip);
        match self.psi_of_psip_interp.as_ref() {
            Some(interp) => Ok(debug_assert_is_finite!(interp.eval_deriv(
                self.psip.uvalues(),
                self.psi.uvalues(),
                psip,
                acc
            )?)),
            None => Err(EvalError::UndefinedEvaluation("dψ(ψp)/dψp".into())),
        }
    }
}

/// Getters.
impl NcQfactor {
    netcdf_path_getter_impl!();
    netcdf_version_getter_impl!();
    equilibrium_type_getter_impl!();
    interp_type_getter_impl!(1);

    /// Returns the value of q on the magnetic axis.
    #[must_use]
    pub fn qaxis(&self) -> f64 {
        match self.q_values.first().copied() {
            Some(qaxis) => qaxis,
            None => unreachable!("NcQfactor cannot be created if `q_values` dont exist"),
        }
    }

    /// Returns the value of q at the last closed flux surface.
    #[must_use]
    pub fn qlast(&self) -> f64 {
        match self.q_values.last().copied() {
            Some(qlast) => qlast,
            None => unreachable!("NcQfactor cannot be created if `q_values` dont exist"),
        }
    }

    lcfs_getter_impl!();
    fluxes_state_getter_impl!();
    fluxes_values_array_getter_impl!();
    vec_to_array1D_getter_impl!(q_array, q_values, q);
}

impl std::fmt::Debug for NcQfactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NcQfactor")
            .field("netCDF path", &self.path())
            .field("netCDF version", &self.netcdf_version().to_string())
            .field("equilibrium type", &self.equilibrium_type())
            .field("interpolation type", &self.interp_type())
            .field("psi", &self.psi)
            .field("psip", &self.psip)
            .field("qaxis", &self.qaxis())
            .field("qlast", &self.qlast())
            .finish()
    }
}

#[cfg(test)]
mod test_utils {
    use super::*;
    use approx::assert_relative_eq;

    pub(super) fn create_nc_qfactor(path_str: &str) -> NcQfactor {
        let path = PathBuf::from(&path_str);
        let builder = NcQfactorBuilder::new(&path, "steffen");
        builder.build().unwrap()
    }

    /// Make sure that dψ(ψp)/dψp and q(ψ) are close enough.
    pub(super) fn test_dpsi_dpsip_q_closeness(qfactor: &impl Qfactor, psip_last: f64, qlast: f64) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psips = Array1::linspace(0.02 * psip_last, 0.98 * psip_last, 100);

        let mut acc = Accelerator::new();
        for psip in psips.iter().copied() {
            assert_relative_eq!(
                qfactor.q_of_psip(psip, &mut acc).unwrap(),
                qfactor.dpsi_dpsip(psip, &mut acc).unwrap(),
                epsilon = qlast * 1e-4
            )
        }
    }

    /// Make sure that dψp(ψ)/dψ and i(ψ) are close enough.
    pub(super) fn test_dpsip_dpsi_iota_closeness(
        qfactor: &impl Qfactor,
        psi_last: f64,
        qlast: f64,
    ) {
        // Do not go to close to the edges, since the interpolation might deviate a bit
        let psis = Array1::linspace(0.02 * psi_last, 0.98 * psi_last, 100);

        let mut acc = Accelerator::new();
        for psi in psis.iter().copied() {
            assert_relative_eq!(
                qfactor.iota_of_psi(psi, &mut acc).unwrap(),
                qfactor.dpsip_dpsi(psi, &mut acc).unwrap(),
                epsilon = qlast * 1e-4
            )
        }
    }
}

#[cfg(test)]
mod test_parabolic_qfactor {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psi_last() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        assert_abs_diff_eq!(qfactor.psi_last(), 0.45);
        assert_relative_eq!(qfactor.psip_last(), 0.25921022097041035, epsilon = 1e-12);
    }

    #[test]
    /// Values calculated at a point where the ParabolicQfactor was defined correctly.
    fn instantiation_with_psip_last() {
        let qfactor = ParabolicQfactor::new(
            1.1,
            3.9,
            LastClosedFluxSurface::Poloidal(0.25921022097041035),
        );
        assert_abs_diff_eq!(qfactor.psip_last(), 0.25921022097041035);
        assert_relative_eq!(qfactor.psi_last(), 0.45, epsilon = 1e-12);
    }
}

#[cfg(test)]
mod test_derivatives_closeness {
    use crate::extract::TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn unity_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = UnityQfactor::new();
        test_dpsi_dpsip_q_closeness(&qfactor, 0.45, 3.9);
    }

    #[test]
    fn unity_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = UnityQfactor::new();
        test_dpsip_dpsi_iota_closeness(&qfactor, 0.45, 3.9);
    }

    #[test]
    fn parabolic_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        let psip_last = qfactor.psip_last();
        let qlast = qfactor.qlast();
        test_dpsi_dpsip_q_closeness(&qfactor, psip_last, qlast);
    }

    #[test]
    fn parabolic_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = ParabolicQfactor::new(1.1, 3.9, LastClosedFluxSurface::Toroidal(0.45));
        let psip_last = qfactor.psip_last();
        let qlast = qfactor.qlast();
        test_dpsip_dpsi_iota_closeness(&qfactor, psip_last, qlast);
    }

    #[test]
    fn nc_qfactor_dpsi_dpsip_q_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        let psip_last = qfactor.psip_last().unwrap();
        let qlast = qfactor.qlast();
        test_dpsi_dpsip_q_closeness(&qfactor, psip_last, qlast);
    }

    #[test]
    fn nc_qfactor_dpsip_dpsi_iota_closeness() {
        let qfactor = create_nc_qfactor(TEST_NETCDF_PATH);
        let psi_last = qfactor.psi_last().unwrap();
        let qlast = qfactor.qlast();
        test_dpsip_dpsi_iota_closeness(&qfactor, psi_last, qlast);
    }
}

#[cfg(test)]
mod test_toroidal_nc_evals {
    use crate::extract::TOROIDAL_TEST_NETCDF_PATH;

    use super::test_utils::*;
    use super::*;

    #[test]
    fn flux_and_interp_states() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        assert_eq!(qfactor.psi_state(), NcFluxState::Good);
        assert_eq!(qfactor.psip_state(), NcFluxState::Bad);
        assert!(qfactor.q_of_psi_interp.is_some());
        assert!(qfactor.q_of_psip_interp.is_none());
        assert!(qfactor.psip_of_psi_interp.is_some());
        assert!(qfactor.psi_of_psip_interp.is_none());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psi_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        assert!(qfactor.q_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.iota_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.psip_of_psi(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.dpsip_dpsi(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psip_evals() {
        let qfactor = create_nc_qfactor(TOROIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psi_of_psip(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsi_dpsip(0.01, &mut acc), Err(err(..)));
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
        assert_eq!(qfactor.psi_state(), NcFluxState::Bad);
        assert_eq!(qfactor.psip_state(), NcFluxState::Good);
        assert!(qfactor.q_of_psi_interp.is_none());
        assert!(qfactor.q_of_psip_interp.is_some());
        assert!(qfactor.psip_of_psi_interp.is_none());
        assert!(qfactor.psi_of_psip_interp.is_some());

        assert!(qfactor.psi_array().is_some());
        assert!(qfactor.psip_array().is_some());
    }

    #[test]
    fn good_psip_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        assert!(qfactor.q_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.iota_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.psi_of_psip(0.01, &mut acc).unwrap().is_finite());
        assert!(qfactor.dpsi_dpsip(0.01, &mut acc).unwrap().is_finite());
    }

    #[test]
    fn bad_psi_evals() {
        let qfactor = create_nc_qfactor(POLOIDAL_TEST_NETCDF_PATH);
        let mut acc = Accelerator::new();
        use EvalError::UndefinedEvaluation as err;
        matches!(qfactor.q_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.iota_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.psip_of_psi(0.01, &mut acc), Err(err(..)));
        matches!(qfactor.dpsip_dpsi(0.01, &mut acc), Err(err(..)));
    }
}
