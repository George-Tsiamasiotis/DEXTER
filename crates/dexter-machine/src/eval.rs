//! Definitions of evaluation methods of machine objects.
//!
//! For analytical machines, this is achieved by evaluation of analytical formulas, while for
//! numerical machines by interpolation over the reconstructed data arrays.

use std::fmt::Debug;

use ndarray::Array1;
use rsl_interpolation::{Accelerator, Accelerator2d};

use crate::{EvalError, FluxCoordinateState, MachineType, MagneticFlux};

/// Reference to a dynamically dispatched [`Mode`] object.
pub type DynMode = Box<dyn Mode>;

/// Reference to a dynamically dispatched [`ModeCache`] object.
pub type DynModeCache = Box<dyn ModeCache + Send + Sync + 'static>;

/// Basic machine object methods.
pub trait MachineObject: Debug + Send + Sync {
    /// Returns the machine object's type.
    fn machine_type(&self) -> MachineType;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψt` flux coordinate.
    fn psit_state(&self) -> FluxCoordinateState;

    /// Returns the [`FluxCoordinateState`] of the toroidal `ψp` flux coordinate.
    fn psip_state(&self) -> FluxCoordinateState;
}

/// Machine geometry related quantities computation.
pub trait Geometry: MachineObject + Debug + Send + Sync {
    /// Returns the magnetic field strength on the axis `B0` in **\[T\]**.
    fn baxis(&self) -> f64;

    /// Returns the horizontal position of the magnetic axis `R0` in **\[m\]**.
    fn raxis(&self) -> f64;

    /// Returns the vertical position of the magnetic axis in **\[m\]**.
    fn zaxis(&self) -> f64;

    /// Returns the horizontal positinon of the geometrical axis (device major radius) **in \[m\]**.
    fn rgeo(&self) -> f64;

    /// Returns the `r` coordinate's value at the last closed flux surface **in \[m\]**.
    fn rlast(&self) -> f64;

    /// Calculates the radial coordinate `r(ψ)` in **\[m\]**.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let r_of_psit: f64 = geometry.eval_r(psit, acc)?;
    /// let r_of_psip: f64 = geometry.eval_r(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_r(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the toroidal flux `ψt(r)` from the radial coordinate `r` **in \[m\]**.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit_of_r: MagneticFlux = geometry.eval_psit_of_r(0.02, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_psit_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError>;

    /// Calculates the poloidal flux `ψp(r)` from the radial coordinate `r` **in \[m\]**.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psip_of_r: MagneticFlux = geometry.eval_psip_of_r(0.02, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_psip_of_r(&self, r: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError>;

    /// Calculates `R(ψ, θ)` **in \[m\]**.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psip = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let rlab_of_psip = geometry.eval_rlab(psip, 3.14, acc)?;
    /// let rlab_of_psip = geometry.eval_rlab(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_rlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `Z(ψ, θ)` **in \[m\]**.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psip = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let zlab_of_psip = geometry.eval_zlab(psip, 3.14, acc)?;
    /// let zlab_of_psip = geometry.eval_zlab(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_zlab(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `J(ψ, θ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let interp1d_type = Interpolation1dType::Akima;
    /// # let interp2d_type = Interpolation2dType::Bicubic;
    /// # let geometry = NcGeometryBuilder::new(&path, interp1d_type, interp2d_type).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let j_of_psit = geometry.eval_jacobian(psit, 3.14, acc)?;
    /// let j_of_psip = geometry.eval_jacobian(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_jacobian(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Returns the last `Rlab` values that correspond to the device's last closed flux surface.
    fn rlab_last(&self) -> Array1<f64>;

    /// Returns the last `Zlab` values that correspond to the device's last closed flux surface.
    fn zlab_last(&self) -> Array1<f64>;
}

/// q-factor related quantities computation.
pub trait Qfactor: MachineObject + Debug + Send + Sync {
    /// Returns the value of the last closed toroidal flux `ψt_last`.
    fn psit_last(&self) -> MagneticFlux;

    /// Returns the value of the last closed poloidal flux `ψp_last`.
    fn psip_last(&self) -> MagneticFlux;

    /// Returns the qfactor's value at the last closed flux surface.
    fn qlast(&self) -> f64;

    /// Returns the qfactor's value at the magnetic axis.
    fn qaxis(&self) -> f64;

    /// Calculates the safety factor `q(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let q_of_psit = qfactor.eval_q(psit, acc)?;
    /// let q_of_psip = qfactor.eval_q(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_q(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the safety factor's derivative `dq(ψ)/dψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let dq_dpsit = qfactor.eval_q_deriv(psit, acc)?;
    /// let dq_dpsip = qfactor.eval_q_deriv(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_q_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Converts a [`MagneticFlux`] to one of the other variant.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the conversion fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Toroidal(0.015);
    ///
    /// let psip_of_psit: MagneticFlux = qfactor.eval_other(psit, acc)?;
    /// let psit_of_psip: MagneticFlux = qfactor.eval_other(psip, acc)?;
    ///
    /// assert_eq!(psit.to_other(&qfactor)?, psip_of_psit);
    /// assert_eq!(psip.to_other(&qfactor)?, psit_of_psip);
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<MagneticFlux, EvalError>;

    /// Calculates `ψt(q)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit_of_q: MagneticFlux = qfactor.eval_psit_of_q(1.2, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_psit_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError>;

    /// Calculates `ψp(q)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psip_of_q: MagneticFlux = qfactor.eval_psip_of_q(1.2, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_psip_of_q(&self, q: f64, acc: &mut Accelerator) -> Result<MagneticFlux, EvalError>;

    /// Calculates the derivative of the **other** [`MagneticFlux`] with respect to `flux`.
    ///
    /// + If a [`MagneticFlux::Toroidal`] is passed, then `dψp/dψt` is calculated.
    /// + If a [`MagneticFlux::Poloidal`] is passed, then `dψt/dψp` is calculated.
    ///
    /// In contrast to [`Qfactor::eval_deriv_wrt_other`], this method only requires one of the
    /// fluxes to be in a "good" state (the one corresponding to the passed `flux` argument).
    ///
    /// This method is useful for ensuring that `dψt/dψp = q` and `dψp/dψt = ι`. The corresponding
    /// methods [`Qfactor::eval_q`] and [`Qfactor::eval_iota`] should be used in calculations as
    /// they are faster and more accurate.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the conversion fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// # use approx::assert_relative_eq;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.1);
    /// let psip = MagneticFlux::Poloidal(0.15);
    ///
    /// let dpsip_dpsit: f64 = qfactor.eval_deriv_of_other(psit, acc)?;
    /// let dpsit_dpsip: f64 = qfactor.eval_deriv_of_other(psip, acc)?;
    ///
    /// let q = qfactor.eval_q(psip, acc)?;
    /// let i = qfactor.eval_iota(psit, acc)?;
    /// assert_relative_eq!(dpsit_dpsip, q);
    /// assert_relative_eq!(dpsip_dpsit, i);
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_of_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<f64, EvalError>;

    /// Calculates a [`MagneticFlux`]'s derivative with respect to the other.
    ///
    /// + If a [`MagneticFlux::Toroidal`] is passed, then `dψt/dψp` is calculated.
    /// + If a [`MagneticFlux::Poloidal`] is passed, then `dψp/dψt` is calculated.
    ///
    /// This method requires **both** fluxes to be in a "good" state. If this is not true,
    /// [`Qfactor::eval_deriv_of_other`] should be used.
    ///
    /// This method is useful for ensuring that `dψt/dψp = q` and `dψp/dψt = ι`. The corresponding
    /// methods [`Qfactor::eval_q`] and [`Qfactor::eval_iota`] should be used in calculations as
    /// they are faster and more accurate.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the conversion fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// # use approx::assert_relative_eq;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.1);
    /// let psip = MagneticFlux::Poloidal(0.15);
    ///
    /// let dpsit_dpsip: f64 = qfactor.eval_deriv_wrt_other(psit, acc)?;
    /// let dpsip_dpsit: f64 = qfactor.eval_deriv_wrt_other(psip, acc)?;
    ///
    /// let q = qfactor.eval_q(psit, acc)?;
    /// let i = qfactor.eval_iota(psip, acc)?;
    /// assert_relative_eq!(dpsit_dpsip, q);
    /// assert_relative_eq!(dpsip_dpsit, i);
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_wrt_other(
        &self,
        flux: MagneticFlux,
        acc: &mut Accelerator,
    ) -> Result<f64, EvalError>;

    /// Calculates the rotational number `ι(ψ) = 1/q(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let iota_of_psit = qfactor.eval_iota(psit, acc)?;
    /// let iota_of_psip = qfactor.eval_iota(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_iota(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        Ok(self.eval_q(flux, acc)?.recip())
    }

    /// Calculates the rotational number's derivative `dι(ψ)/dψ = -q'/q²`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let lcfs = MagneticFlux::Toroidal(0.45);
    /// # let qfactor = ParabolicQfactor::new(1.1, 3.8, lcfs);
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let diota_dpsit = qfactor.eval_iota_deriv(psit, acc)?;
    /// let diota_dpsip = qfactor.eval_iota_deriv(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_iota_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError> {
        Ok(-self.eval_q_deriv(flux, acc)? / self.eval_q(flux, acc)?.powi(2))
    }
}

/// Plasma current related quantities computation.
pub trait Current: MachineObject + Debug + Send + Sync {
    /// Calculates the poloidal current `g(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let g_of_psit = current.eval_g(psit, acc)?;
    /// let g_of_psip = current.eval_g(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_g(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the toroidal current `I(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let i_of_psit = current.eval_i(psit, acc)?;
    /// let i_of_psip = current.eval_i(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_i(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the poloidal current's derivative `dg/dψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let dg_dpsit = current.eval_g_deriv(psit, acc)?;
    /// let dg_dpsip = current.eval_g_deriv(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_g_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the toroidal current's derivative `dI/dψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let di_dpsit = current.eval_i_deriv(psit, acc)?;
    /// let di_dpsip = current.eval_i_deriv(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_i_deriv(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the poloidal current's second derivative `d²g/dψ²`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let dg_dpsit2 = current.eval_g_deriv2(psit, acc)?;
    /// let dg_dpsip2 = current.eval_g_deriv2(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_g_deriv2(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;

    /// Calculates the toroidal current's second derivative `d²I/dψ²`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use std::path::PathBuf;
    /// # use rsl_interpolation::Accelerator;
    /// #
    /// # let current = LarCurrent;
    /// #
    /// let acc = &mut Accelerator::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let di_dpsit2 = current.eval_i_deriv2(psit, acc)?;
    /// let di_dpsip2 = current.eval_i_deriv2(psip, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_i_deriv2(&self, flux: MagneticFlux, acc: &mut Accelerator) -> Result<f64, EvalError>;
}

/// Magnetic field related quantities computation.
pub trait Bfield: MachineObject + Debug + Send + Sync {
    /// Calculates `B(ψ, θ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let b_of_psit = bfield.eval_b(psit, 3.14, acc)?;
    /// let b_of_psip = bfield.eval_b(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_b(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψ, θ)/dψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let db_dpsit = bfield.eval_deriv_flux(psit, 3.14, acc)?;
    /// let db_dpsip = bfield.eval_deriv_flux(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_flux(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `dB(ψ, θ)/dθ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let db_of_psit_dtheta = bfield.eval_deriv_theta(psit, 3.14, acc)?;
    /// let db_of_psip_dtheta = bfield.eval_deriv_theta(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_theta(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `d²B(ψ, θ)/dψ²`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let db_dpsit2 = bfield.eval_deriv_flux2(psit, 3.14, acc)?;
    /// let db_dpsip2 = bfield.eval_deriv_flux2(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_flux2(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `d²B(ψ, θ)/dθ²`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let db_of_psit_dtheta2 = bfield.eval_deriv_theta2(psit, 3.14, acc)?;
    /// let db_of_psip_dtheta2 = bfield.eval_deriv_theta2(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_theta2(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;

    /// Calculates `d²B(ψ, θ)/dθdψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # use std::path::PathBuf;
    /// #
    /// # let path = PathBuf::from("./netcdf.nc");
    /// # let bfield = NcBfieldBuilder::new(&path, Interpolation2dType::Bilinear).build()?;
    /// #
    /// let acc = &mut Accelerator2d::new();
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let psip = MagneticFlux::Poloidal(0.015);
    ///
    /// let db_dpsit_dtheta = bfield.eval_deriv_mixed(psit, 3.14, acc)?;
    /// let db_dpsip_dtheta = bfield.eval_deriv_mixed(psip, 3.14, acc)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_mixed(
        &self,
        flux: MagneticFlux,
        theta: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, EvalError>;
}

/// Electric field related quantities computation.
pub trait Efield: MachineObject + Debug + Send + Sync {
    /// Calculates `E(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let efield = GaussianEfield::new(5e-4, lcfs * 0.9, 0.02 * lcfs.value());
    ///
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let e = efield.eval_e(psit)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_e(&self, flux: MagneticFlux) -> Result<f64, EvalError>;

    /// Calculates `Φ(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Poloidal(0.35);
    /// let efield = GaussianEfield::new(5e-4, lcfs * 0.9, 0.02 * lcfs.value());
    ///
    /// let psip = MagneticFlux::Poloidal(0.01);
    /// let phi = efield.eval_phi(psip)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_phi(&self, flux: MagneticFlux) -> Result<f64, EvalError>;

    /// Calculates `dΦ(ψ)/dψ`.
    ///
    /// Since the field is radial, this method is equivalent to `-E(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let efield = GaussianEfield::new(5e-4, lcfs * 0.9, 0.02 * lcfs.value());
    ///
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let dphi = efield.eval_phi_deriv(psit)?;
    ///
    /// assert_eq!(dphi, -efield.eval_e(psit)?);
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_phi_deriv(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        Ok(-self.eval_e(flux)?)
    }

    /// Calculates `d²Φ(ψ)/dψ²`.
    ///
    /// Since the field is radial, this method is equivalent to `-dE(ψ)/d(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let efield = GaussianEfield::new(5e-4, lcfs * 0.9, 0.02 * lcfs.value());
    ///
    /// let psit = MagneticFlux::Toroidal(0.01);
    /// let dphi2 = efield.eval_phi_deriv2(psit)?;
    ///
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_phi_deriv2(&self, flux: MagneticFlux) -> Result<f64, EvalError>;
}

/// Defines the behavior of objects that support caching of a [`Mode`]'s values.
///
/// All non-trivial [`Mode`] evaluation methods must be guarded by a:
/// ```notest
/// if !cache.is_updated() {
///     cache.update()
/// }
/// ```
/// statement before using the cache.
#[expect(
    private_bounds,
    reason = "only used internally for creating Perturbation"
)]
pub trait ModeCache: DynModeCacheClone + Debug {
    /// Checks if the cache's stored independent coordinates are up-to-date, i.e. are equal to the
    /// passed arguments.
    fn is_updated(&mut self, flux: f64, theta: f64, zeta: f64, t: f64) -> bool;

    /// Updates the cache's coordinates and intermediate values.
    fn update(&mut self, flux: f64, theta: f64, zeta: f64, t: f64);

    /// Returns a mutable reference to the cache's `cache` array.
    fn cache(&mut self) -> &mut [f64];

    /// Returns a reference to the cache's `params` array.
    fn params(&mut self) -> &[f64];

    /// Returns a mutable reference to the cache's [`Accelerator`], if it exists.
    fn acc(&mut self) -> Option<&mut Accelerator>;

    /// Returns the cache's hits.
    fn hits(&self) -> usize;

    /// Returns the cache's misses.
    fn misses(&self) -> usize;
}

/// Single perturbation mode related quantities computation.
#[expect(
    private_bounds,
    reason = "only used internally for creating Perturbation"
)]
pub trait Mode: MachineObject + DynModeClone + Debug + Send + Sync {
    /// Returns the poloidal mode number `m`.
    fn m(&self) -> i64;

    /// Returns the toroidal mode number `n`.
    fn n(&self) -> i64;

    /// Returns a default instance of the Mode's corresponding caching object.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    /// let cache = mode.generate_cache();
    /// # Ok::<_, MachineError>(())
    /// ```
    fn generate_cache(&self) -> DynModeCache;

    /// Calculates the mode's amplitude `α(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psit = MagneticFlux::Toroidal(0.01);
    ///
    /// let a_of_psit = mode.eval_amplitude(psit, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_amplitude(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's phase `φ(ψ)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let phase_of_psip = mode.eval_phase(psip, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_phase(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's value `m(ψ, θ, ζ, t)`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let m_of_psip = mode.eval_m(psip, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_m(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to the magnetic flux `dm(ψ, θ, ζ, t)/dψ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psit = MagneticFlux::Toroidal(0.01);
    ///
    /// let dm_dpsit = mode.eval_deriv_flux(psit, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_flux(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to `θ`, `dm(ψ, θ, ζ, t)/dθ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psit = MagneticFlux::Toroidal(0.01);
    ///
    /// let dm_dtheta = mode.eval_deriv_theta(psit, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_theta(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to `ζ`, `dm(ψ, θ, ζ, t)/dζ`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let dm_dzeta = mode.eval_deriv_zeta(psip, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_zeta(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;

    /// Calculates the mode's derivative with respect to `t`, `dm(ψ, θ, ζ, t)/dt`.
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the evaluation fails for any reason.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Poloidal(0.45);
    /// let mode = FluteMode::new(1e-3, lcfs, 3, 2, 0.0);
    ///
    /// let mut cache = mode.generate_cache();
    /// let psip = MagneticFlux::Poloidal(0.01);
    ///
    /// let dm_dt = mode.eval_deriv_t(psip, 0.2, 0.3, 0.0, &mut cache)?;
    /// # Ok::<_, MachineError>(())
    /// ```
    fn eval_deriv_t(
        &self,
        flux: MagneticFlux,
        theta: f64,
        zeta: f64,
        t: f64,
        cache: &mut DynModeCache,
    ) -> Result<f64, EvalError>;
}

// ================================================================================================

// HACK: This is necessary to clone `dynMode` = `Box<dyn Mode>`.
// https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object

/// Clone `Box<dyn Mode>`.
trait DynModeClone {
    /// Clone box.
    fn clone_box(&self) -> DynMode;
}

/// Clone `Box<dyn ModeCache>`.
trait DynModeCacheClone {
    /// Clone box.
    fn clone_box(&self) -> DynModeCache;
}

impl<M> DynModeClone for M
where
    M: 'static + Mode + Clone,
{
    fn clone_box(&self) -> DynMode {
        Box::new(self.clone())
    }
}

impl<C> DynModeCacheClone for C
where
    C: 'static + ModeCache + Clone + Send + Sync + 'static,
{
    fn clone_box(&self) -> DynModeCache {
        Box::new(self.clone())
    }
}

impl Clone for DynMode {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

impl Clone for DynModeCache {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
