//! Defines the `MagneticFluxKind` and `MagneticFlux` enums.

use ndarray::Array1;
use rsl_interpolation::Accelerator;

use crate::{EvalError, MagneticFluxArray, Qfactor};

/// The kind of a magnetic flux coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MagneticFluxKind {
    /// Toroidal magnetic flux `ψt`.
    Toroidal,
    /// Poloidal magnetic flux `ψp`.
    Poloidal,
}

impl MagneticFluxKind {
    /// Creates a [`MagneticFlux`] of `self`'s kind and value `value`.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// assert_eq!(
    ///     MagneticFluxKind::Toroidal.to_magnetic_flux(0.1),
    ///     MagneticFlux::Toroidal(0.1),
    /// );
    /// ```
    #[must_use]
    pub fn to_magnetic_flux(&self, value: f64) -> MagneticFlux {
        match *self {
            Self::Toroidal => MagneticFlux::Toroidal(value),
            Self::Poloidal => MagneticFlux::Poloidal(value),
        }
    }

    /// Creates a [`MagneticFluxArray`] of `self`'s kind and values `values`.
    #[must_use]
    pub fn to_magnetic_flux_array(&self, values: Array1<f64>) -> MagneticFluxArray {
        match *self {
            Self::Toroidal => MagneticFluxArray::Toroidal(values),
            Self::Poloidal => MagneticFluxArray::Poloidal(values),
        }
    }

    /// Returns the symbol of `self` ("ψt" or "ψp") as a `Box<str>`.
    #[must_use]
    pub fn symbol(&self) -> Box<str> {
        match *self {
            Self::Toroidal => "ψt".into(),
            Self::Poloidal => "ψp".into(),
        }
    }
}

/// Representation of the magnetic flux and its kind.
#[derive(Clone, Copy, PartialEq)]
pub enum MagneticFlux {
    /// Toroidal magnetic flux `ψt`.
    Toroidal(f64),
    /// Poloidal magnetic flux `ψp`.
    Poloidal(f64),
}

impl MagneticFlux {
    /// Returns the value of `self`, regardless of the kind.
    #[must_use]
    pub fn value(&self) -> f64 {
        match *self {
            Self::Toroidal(value) | Self::Poloidal(value) => value,
        }
    }

    /// Returns a mutable reference to the value of `self`, regardless the kind.
    #[must_use]
    pub fn value_mut(&mut self) -> &mut f64 {
        match *self {
            Self::Toroidal(ref mut value) | Self::Poloidal(ref mut value) => value,
        }
    }

    /// Returns the symbol of `self` ("ψt" or "ψp") as a `Box<str>`.
    #[must_use]
    pub fn symbol(&self) -> Box<str> {
        match *self {
            Self::Toroidal(_) => "ψt".into(),
            Self::Poloidal(_) => "ψp".into(),
        }
    }

    /// Returns the kind of `self`.
    #[must_use]
    pub fn kind(&self) -> MagneticFluxKind {
        match *self {
            Self::Toroidal(_) => MagneticFluxKind::Toroidal,
            Self::Poloidal(_) => MagneticFluxKind::Poloidal,
        }
    }

    /// Creates a new `MagneticFlux` where the inner value has been converted to the other magnetic
    /// flux through [`Qfactor::eval_other`][crate::Qfactor::eval_other].
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the contained value is out-of-bounds or the conversion is not defined.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # use rsl_interpolation::*;
    /// # let acc = &mut Accelerator::new();
    /// let qfactor = ParabolicQfactor::new(1.1, 3.9, MagneticFlux::Toroidal(0.5));
    /// let psit = MagneticFlux::Toroidal(0.1);
    /// let psip = psit.to_other(&qfactor)?;
    /// assert_eq!(psip.kind(), MagneticFluxKind::Poloidal);
    /// assert_eq!(psip, qfactor.eval_other(psit, acc)?);
    /// # Ok::<_, EvalError>(())
    /// ```
    pub fn to_other(&self, qfactor: &dyn Qfactor) -> Result<Self, EvalError> {
        qfactor.eval_other(*self, &mut Accelerator::new())
    }

    /// Converts the magnetic flux into the other kind through
    /// [`Qfactor::eval_other`][crate::Qfactor::eval_other].
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if the contained value is out-of-bounds or the conversion is not defined.
    pub fn into_other(mut self, qfactor: &dyn Qfactor) -> Result<Self, EvalError> {
        self = self.to_other(qfactor)?;
        Ok(self)
    }
}

impl std::fmt::Debug for MagneticFlux {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.symbol(), self.value())
    }
}

// ===============================================================================================

impl std::ops::Add<f64> for MagneticFlux {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        use MagneticFlux::*;
        match self {
            Toroidal(value) => Toroidal(value + rhs),
            Poloidal(value) => Poloidal(value + rhs),
        }
    }
}

impl std::ops::Sub<f64> for MagneticFlux {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        use MagneticFlux::*;
        match self {
            Toroidal(value) => Toroidal(value - rhs),
            Poloidal(value) => Poloidal(value - rhs),
        }
    }
}

impl std::ops::Mul<f64> for MagneticFlux {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        use MagneticFlux::*;
        match self {
            Toroidal(value) => Toroidal(value * rhs),
            Poloidal(value) => Poloidal(value * rhs),
        }
    }
}

impl std::ops::Div<f64> for MagneticFlux {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        use MagneticFlux::*;
        match self {
            Toroidal(value) => Toroidal(value / rhs),
            Poloidal(value) => Poloidal(value / rhs),
        }
    }
}

impl std::ops::Neg for MagneticFlux {
    type Output = Self;

    fn neg(self) -> Self::Output {
        use MagneticFlux::*;
        match self {
            Toroidal(value) => Toroidal(-value),
            Poloidal(value) => Poloidal(-value),
        }
    }
}

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::AbsDiffEq for MagneticFlux {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::EPSILON
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        (self.value() - other.value()).abs() <= epsilon
    }
}

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::RelativeEq for MagneticFlux {
    fn default_max_relative() -> Self::Epsilon {
        f64::EPSILON
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        self.value()
            .relative_eq(&other.value(), epsilon, max_relative)
    }
}

/// Compares discriminants and then forwards the comparison to the contained `f64`.
impl approx::UlpsEq for MagneticFlux {
    fn default_max_ulps() -> u32 {
        4
    }

    fn ulps_eq(&self, other: &Self, epsilon: Self::Epsilon, max_ulps: u32) -> bool {
        use std::mem::discriminant;
        if discriminant(self) != discriminant(other) {
            return false;
        }
        self.value().ulps_eq(&other.value(), epsilon, max_ulps)
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::{MagneticFlux::*, ParabolicQfactor};
    use approx::*;
    use ndarray::array;

    #[test]
    fn magnetic_flux_kind_to_magnetic_flux() {
        assert_eq!(
            MagneticFluxKind::Toroidal.to_magnetic_flux(0.1),
            Toroidal(0.1)
        );
        assert_eq!(
            MagneticFluxKind::Poloidal.to_magnetic_flux(0.2),
            Poloidal(0.2)
        );
    }

    #[test]
    fn magnetic_flux_kind_to_magnetic_flux_array() {
        let f = MagneticFluxKind::Toroidal.to_magnetic_flux_array(array![0.1, 0.2]);
        assert_eq!(f.kind(), MagneticFluxKind::Toroidal);
        assert_eq!(array![0.1, 0.2], f.values());

        let f = MagneticFluxKind::Poloidal.to_magnetic_flux_array(array![0.3, 0.4]);
        assert_eq!(f.kind(), MagneticFluxKind::Poloidal);
        assert_eq!(array![0.3, 0.4], f.values());
    }

    #[test]
    fn add() {
        assert_relative_eq!(Toroidal(0.1) + 0.4, Toroidal(0.5));
        assert_relative_eq!(Poloidal(0.1) + 0.5, Poloidal(0.6));
    }

    #[test]
    fn sun() {
        assert_relative_eq!(Toroidal(0.4) - 0.1, Toroidal(0.3));
        assert_relative_eq!(Poloidal(0.5) - 0.1, Poloidal(0.4));
    }

    #[test]
    fn mul() {
        assert_relative_eq!(Toroidal(0.1) * 0.4, Toroidal(0.04));
        assert_relative_eq!(Poloidal(0.1) * 0.5, Poloidal(0.05));
    }

    #[test]
    fn div() {
        assert_relative_eq!(Toroidal(0.4) / 0.2, Toroidal(2.0));
        assert_relative_eq!(Poloidal(0.4) / 0.2, Poloidal(2.0));
    }

    #[test]
    fn neg() {
        assert_relative_eq!(-Toroidal(0.4), Toroidal(-0.4));
        assert_relative_eq!(-Poloidal(0.4), Poloidal(-0.4));
    }

    #[test]
    fn relative_eq() {
        assert_relative_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_relative_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_relative_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_relative_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_relative_ne!(Toroidal(0.01), Poloidal(0.01));
    }

    #[test]
    fn abs_diff_eq() {
        assert_abs_diff_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_abs_diff_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_abs_diff_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_abs_diff_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_abs_diff_ne!(Toroidal(0.01), Poloidal(0.01));
    }

    #[test]
    fn ulps_eq() {
        assert_ulps_eq!(Toroidal(0.01), Toroidal(0.01));
        assert_ulps_eq!(Poloidal(0.01), Poloidal(0.01));

        assert_ulps_ne!(Toroidal(0.01), Toroidal(0.02));
        assert_ulps_ne!(Poloidal(0.01), Poloidal(0.02));
        assert_ulps_ne!(Toroidal(0.01), Poloidal(0.01));
    }

    #[test]
    fn to_other() {
        let lcfs = Toroidal(1.0);
        let qfactor = ParabolicQfactor::new(1.0, 2.0, lcfs); // q = ψt^2 + 1

        let tor = Toroidal(0.3);
        let other_pol = tor.to_other(&qfactor).unwrap();

        let pol = Poloidal(other_pol.value());
        let other_tor = pol.to_other(&qfactor).unwrap();

        assert_relative_eq!(tor, other_tor);
        assert_relative_eq!(pol, other_pol);
    }
}
