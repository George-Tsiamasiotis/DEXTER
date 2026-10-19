//! Defines the `MagneticFluxArray` enum.

use ndarray::{Array1, ArrayView1, ArrayViewMut1};
use rsl_interpolation::Accelerator;

use crate::{EvalError, MagneticFlux, MagneticFluxKind, Qfactor};

/// Representation of an [`Array1`] of magnetic flux values and their kind.
#[derive(Clone, PartialEq)]
pub enum MagneticFluxArray {
    /// Toroidal magnetic flux `ψ`.
    Toroidal(Array1<f64>),
    /// Poloidal magnetic flux `ψp`.
    Poloidal(Array1<f64>),
}

impl std::ops::Index<usize> for MagneticFluxArray {
    type Output = f64;
    /// Returns the indexed value of the contained array, regardless of the kind.
    fn index(&self, index: usize) -> &Self::Output {
        match *self {
            Self::Toroidal(ref values) | Self::Poloidal(ref values) => &values[index],
        }
    }
}

impl MagneticFluxArray {
    /// Returns a view on the values of `self`, regardless of the kind.
    #[must_use]
    pub fn values(&self) -> ArrayView1<'_, f64> {
        match *self {
            Self::Toroidal(ref values) | Self::Poloidal(ref values) => values.view(),
        }
    }

    /// Returns a mutable view on the values of `self`, regardless of the kind.
    #[must_use]
    pub fn values_mut(&mut self) -> ArrayViewMut1<'_, f64> {
        match *self {
            Self::Toroidal(ref mut values) | Self::Poloidal(ref mut values) => values.view_mut(),
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

    /// Returns the symbol of `self` ("ψ" or "ψp") as a `Box<str>`.
    #[must_use]
    pub fn symbol(&self) -> Box<str> {
        self.kind().symbol()
    }

    /// Returns the number of points of the inner array.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values().len()
    }

    /// Returns `true` if `self` has no values.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns a [`MagneticFlux`] of kind `self.kind()` and value `self.values()[index]`.
    ///
    /// Returns `None` if `index` is greater than `self.len()`.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<MagneticFlux> {
        if index >= self.len() {
            return None;
        }
        match *self {
            Self::Toroidal(ref values) => Some(MagneticFlux::Toroidal(values[index])),
            Self::Poloidal(ref values) => Some(MagneticFlux::Poloidal(values[index])),
        }
    }

    /// Sets the `index`-th element of the array to `value`, regardless of the kind.
    ///
    /// # Panics
    ///
    /// This method panics if `index` is greater than `self.len()`.
    pub fn set_value(&mut self, index: usize, value: f64) {
        self.values_mut()[index] = value;
    }

    /// Creates a new `MagneticFluxArray` where all values have been converted to the other magnetic
    /// flux through [`Qfactor::eval_other`][crate::Qfactor::eval_other].
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the contained values are out-of-bounds or the conversion
    /// is not defined.
    pub fn to_other(&self, qfactor: &dyn Qfactor) -> Result<Self, EvalError> {
        let acc = &mut Accelerator::new();
        let mut new_values = Vec::<f64>::with_capacity(self.len());

        for index in 0..self.len() {
            let Some(flux) = self.get(index) else {
                unreachable!("index is in-bounds")
            };
            new_values.push(qfactor.eval_other(flux, acc)?.value());
        }

        match self.kind() {
            MagneticFluxKind::Toroidal => Ok(Self::Poloidal(Array1::from_vec(new_values))),
            MagneticFluxKind::Poloidal => Ok(Self::Toroidal(Array1::from_vec(new_values))),
        }
    }

    /// Converts the magnetic flux array into the other kind through
    /// [`Qfactor::eval_other`][crate::Qfactor::eval_other].
    ///
    /// # Errors
    ///
    /// Returns an [`EvalError`] if any of the contained values are out-of-bounds or the conversion
    /// is not defined.
    pub fn into_other(mut self, qfactor: &dyn Qfactor) -> Result<Self, EvalError> {
        self = self.to_other(qfactor)?;
        Ok(self)
    }
}

impl std::fmt::Debug for MagneticFluxArray {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.symbol(), self.values())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ParabolicQfactor;

    #[test]
    fn magnetic_flux_array() {
        let values = Array1::from_vec(vec![0.1, 0.2]);
        let arr = MagneticFluxArray::Toroidal(values.clone());
        assert_eq!(arr.kind(), MagneticFluxKind::Toroidal);
        assert_eq!(arr.values(), values);
        assert_eq!(arr.len(), 2);
        assert!(!arr.is_empty());
        assert_eq!(arr.get(0).unwrap(), MagneticFlux::Toroidal(0.1));
        assert_eq!(arr.get(1).unwrap(), MagneticFlux::Toroidal(0.2));

        let values = Array1::from_vec(vec![0.1, 0.2]);
        let arr = MagneticFluxArray::Poloidal(values.clone());
        assert_eq!(arr.kind(), MagneticFluxKind::Poloidal);
        assert_eq!(arr.values(), values);
        assert_eq!(arr.len(), 2);
        assert!(!arr.is_empty());
        assert_eq!(arr.get(0).unwrap(), MagneticFlux::Poloidal(0.1));
        assert_eq!(arr.get(1).unwrap(), MagneticFlux::Poloidal(0.2));
    }

    #[test]
    fn magnetic_flux_array_to_other() {
        let lcfs = MagneticFlux::Toroidal(1.0);
        let qfactor = ParabolicQfactor::new(1.0, 2.0, lcfs); // q = ψ^2 + 1
        let values = Array1::from_vec(vec![0.1, 0.2]);
        let tor = MagneticFluxArray::Toroidal(values.clone());

        let other_pol = tor.to_other(&qfactor).unwrap();
        assert_eq!(other_pol.kind(), MagneticFluxKind::Poloidal);
        assert_eq!(other_pol.len(), 2);
        assert!(!other_pol.is_empty());
        assert_eq!(other_pol.get(0).unwrap().kind(), MagneticFluxKind::Poloidal);
        assert_eq!(other_pol.get(1).unwrap().kind(), MagneticFluxKind::Poloidal);

        let pol = MagneticFluxArray::Poloidal(other_pol.values().to_owned());

        let other_tor = pol.to_other(&qfactor).unwrap();
        assert_eq!(other_tor.kind(), MagneticFluxKind::Toroidal);
        assert_eq!(other_tor.len(), 2);
        assert!(!other_tor.is_empty());
        assert_eq!(other_tor.get(0).unwrap().kind(), MagneticFluxKind::Toroidal);
        assert_eq!(other_tor.get(1).unwrap().kind(), MagneticFluxKind::Toroidal);

        assert!(tor.values().relative_eq(&other_tor.values(), 1e-10, 1e-10));
        assert!(pol.values().relative_eq(&other_pol.values(), 1e-10, 1e-10));
    }

    #[test]
    fn magnetic_flux_array_set_values() {
        let values = Array1::from_vec(vec![0.1, 0.2]);
        let mut flux = MagneticFluxArray::Toroidal(values.clone());

        flux.set_value(1, 0.3);
        assert_eq!(flux[1], 0.3);
    }
}
