//! Representation of a machine's electric field.

use statrs::function::erf::erf;

use crate::{Efield, EvalError, FluxCoordinateState, MachineObject, MachineType, MagneticFlux};
use crate::{MagneticFluxKind, debug_assert_non_negative_flux};

// ===============================================================================================

/// No electric field.
///
/// This type does not need to be constructed directly; It is created when
/// [`Machine`][`crate::Machine`] is built without calling
/// [`MachineBuilder::with_efield`][crate::MachineBuilder::with_efield].
///
/// # Note
///
/// No ψ/ψp bounds checks are performed in evaluations.
#[derive(Debug)]
pub struct NoEfield;

impl MachineObject for NoEfield {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }

    fn psip_state(&self) -> FluxCoordinateState {
        FluxCoordinateState::Good
    }
}

impl Efield for NoEfield {
    fn eval_e(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }

    fn eval_phi(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }

    fn eval_phi_deriv(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        Ok(0.0)
    }
}

// ===============================================================================================

/// Gaussian electric field of the form `E(ψ) = Ε⋅exp[-(ψ - μ)²/(2σ²)]`, where:
/// + `E` corresponds to the electric field strength at the curve's peak, in **Normalized Units**,
/// + `μ` is the magnetic flux surface where the peak is located,
/// + `σ` is the *standard deviation* or *root mean square width* of the 'bell'.
///
/// # Note
///
/// + The [`kind`][MagneticFluxKind] of `mean` defines the electric field's [`Good`] flux coordinate.
/// + No ψ/ψp bounds checks are performed in evaluations.
///
/// [`Good`]: FluxCoordinateState::Good
#[non_exhaustive]
pub struct GaussianEfield {
    /// The electric field strength at the curve's peak **in Normalized Units**.
    strength: f64,
    /// The magnetic flux surface where the peak is located.
    mean: MagneticFlux,
    /// The *standard deviation* or *root mean square width* of the 'bell'.
    width: f64,
    /// The kind of the flux w.r.t. the field is defined.
    _good_flux_kind: MagneticFluxKind,
    /// (Intermediate value) The variance σ².
    _var: f64,
    /// (Intermediate value) The coefficient `E⋅σ⋅sqrt(π/2)`.
    _strength_coef: f64,
    /// (Intermediate value) The constant `erf(μ/(sqrt(2)⋅σ))`.
    _erf_const_term: f64,
    /// (Intermediate value) The constant `sqrt(2)/(2σ)`.
    _erf_coef: f64,
}

impl GaussianEfield {
    /// Creates a new `GaussianEfield`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    ///
    /// let strength = 5e-4; // Normalized units
    /// let mean = lcfs * 0.9; // Peak close to the LCFS
    /// let width = 0.02 * lcfs.value(); // width of 2% of the ψ span.
    ///
    /// let efield = GaussianEfield::new(strength, mean, width);
    /// ```
    #[must_use]
    pub fn new(strength: f64, mean: MagneticFlux, width: f64) -> Self {
        use std::f64::consts::{FRAC_PI_2, SQRT_2};
        let _good_flux_kind = mean.kind();
        let _var = width.powi(2);
        let _strength_coef = strength * width * FRAC_PI_2.sqrt();
        let _erf_const_term = erf(mean.value() / (SQRT_2 * width));
        let _erf_coef = SQRT_2 / (2.0 * width);

        Self {
            strength,
            mean,
            width,
            _good_flux_kind,
            _var,
            _strength_coef,
            _erf_const_term,
            _erf_coef,
        }
    }
}

impl MachineObject for GaussianEfield {
    fn machine_type(&self) -> MachineType {
        MachineType::Analytical
    }

    fn psi_state(&self) -> FluxCoordinateState {
        match self.mean.kind() {
            MagneticFluxKind::Toroidal => FluxCoordinateState::Good,
            MagneticFluxKind::Poloidal => FluxCoordinateState::Bad,
        }
    }

    fn psip_state(&self) -> FluxCoordinateState {
        match self.mean.kind() {
            MagneticFluxKind::Toroidal => FluxCoordinateState::Bad,
            MagneticFluxKind::Poloidal => FluxCoordinateState::Good,
        }
    }
}

impl Efield for GaussianEfield {
    fn eval_e(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.kind() != self._good_flux_kind {
            let msg = format!("E({})", flux.kind().symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        }
        let mean = self.mean.value();
        let width = self.width;
        let flux = flux.value();

        let exponent = -0.5 * ((flux - mean) / width).powi(2);
        Ok(-self.strength * exponent.exp())
    }

    fn eval_phi(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.kind() != self._good_flux_kind {
            let msg = format!("Φ({})", flux.kind().symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        }
        let mean = self.mean.value();
        let flux = flux.value();

        let erf_arg = self._erf_coef * (flux - mean);
        Ok(self._strength_coef * (erf(erf_arg) + self._erf_const_term))
    }

    fn eval_phi_deriv(&self, flux: MagneticFlux) -> Result<f64, EvalError> {
        debug_assert_non_negative_flux!(flux);
        if flux.kind() != self._good_flux_kind {
            let msg = format!("dΦ({})/d{}", flux.kind().symbol(), flux.kind().symbol());
            return Err(EvalError::UndefinedEvaluation(msg));
        }
        Ok(-self.eval_e(flux)?)
    }
}

/// Getters.
impl GaussianEfield {
    /// Returns the electric field strength at the curve's peak **in Normalized Units**.
    #[must_use]
    pub fn strength(&self) -> f64 {
        self.strength
    }

    /// Returns the magnetic flux surface where the peak is located.
    #[must_use]
    pub fn mean(&self) -> MagneticFlux {
        self.mean
    }

    /// Returns the *standard deviation* or *root mean square width* of the 'bell'.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.width
    }
}

impl std::fmt::Debug for GaussianEfield {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GaussianEfield")
            .field("strength", &self.strength)
            .field("mean", &self.mean)
            .field("width", &self.width)
            .finish()
    }
}
