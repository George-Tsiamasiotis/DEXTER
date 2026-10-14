//! Defines [`ParametrizedSegment`], a container type to store and parametrize the Stationary
//! Curve's distinct segments.

use dexter_machine::{Machine, MagneticFluxArray, MagneticFluxKind};
use ndarray::{Array1, ArrayView1};

use crate::{COMError, GCHamiltonian};

/// A 2D parametrized curve. Represents an individual segment of the
/// [`StationaryCurve`][crate::StationaryCurve].
///
/// The parameter `t` is chosen to be the index of the corresponding `(ψ,θ)` point.
#[derive(Clone)]
pub struct ParametrizedSegment {
    /// The segment's parameter.
    t: Array1<f64>,
    /// The segment's magnetic flux values.
    flux_array: MagneticFluxArray,
    /// The segment's `θ` values.
    theta_array: Array1<f64>,
    /// The segment's `dot()θ)` values, calculated on each `(ψ, θ)` pair of the `flux` and `theta`
    /// arrays.
    theta_dot: Option<Array1<f64>>,
}

impl ParametrizedSegment {
    /// Creates a new `ParametrizedSegment` from the `flux` and `theta` array and initializes the
    /// parameter `t`.
    ///
    /// # Panics
    ///
    /// This method panics if the input arrays do not have the same length.
    #[must_use]
    pub(crate) fn new(flux_array: MagneticFluxArray, theta_array: Array1<f64>) -> Self {
        assert_eq!(
            flux_array.len(),
            theta_array.len(),
            "Input arrays must have the same length"
        );
        let t = Array1::linspace(0.0, (flux_array.len() - 1) as f64, flux_array.len());
        Self {
            t,
            flux_array,
            theta_array,
            theta_dot: None,
        }
    }

    /// Returns the number of points that comprise the segment.
    #[must_use]
    pub fn len(&self) -> usize {
        self.t.len()
    }

    /// Returns `true` if the segment has no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the [`MagneticFluxArray`] of the stored [`MagneticFluxArray`].
    #[must_use]
    pub fn flux_kind(&self) -> MagneticFluxKind {
        self.flux_array.kind()
    }

    /// Returns a view on the `t` parameter array.
    #[must_use]
    pub fn t_array(&self) -> ArrayView1<'_, f64> {
        self.t.view()
    }

    /// Returns a view on the contained [`MagneticFluxArray`].
    #[must_use]
    pub fn flux_array(&self) -> &MagneticFluxArray {
        &self.flux_array
    }

    /// Returns a view on the `θ` values array.
    #[must_use]
    pub fn theta_array(&self) -> ArrayView1<'_, f64> {
        self.theta_array.view()
    }
    /// Returns a view on the `dot(θ)` values array, if they have been calculated.
    #[must_use]
    pub fn theta_dot_array(&self) -> Option<ArrayView1<'_, f64>> {
        self.theta_dot.as_ref().map(|array| array.view())
    }
}

/// Calculations.
impl ParametrizedSegment {
    /// Calculates the `dot(θ)` array from the contained `flux` and `theta` arrays and the constants
    /// of motion `Pζ` and `μ`.
    ///
    /// # Errors
    ///
    /// + Returns a [`COMError::BadFluxCoordinate`] if `self.flux_array.kind()` is not a [`Good`]
    ///   coordinate of `machine`.
    /// + Returns a [`COMError::EvalError`] if `flux_array` contains any out-of-bounds values or
    ///   if any evaluation fails.
    ///
    /// [`Good`]: dexter_machine::FluxCoordinateState::Good
    pub fn calculate_theta_dot(
        &mut self,
        machine: Machine,
        pzeta: f64,
        mu: f64,
    ) -> Result<(), COMError> {
        self.theta_dot = Some(GCHamiltonian::eval_theta_dot_array(
            machine,
            &self.theta_array,
            &self.flux_array,
            pzeta,
            mu,
        )?);
        Ok(())
    }
}

impl std::fmt::Debug for ParametrizedSegment {
    #[expect(clippy::unwrap_used, reason = "safe")]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // NOTE: By construction, `self` always has a length of at least 2, therefore `unwrap` is
        // safe here.
        let last_index = self.len() - 1;
        let flux_last = self.flux_array.get(last_index).unwrap().value();
        let start = format_args!(
            "flux(normalized)={:?}, θ={}",
            self.flux_array.get(0).unwrap() / flux_last,
            self.theta_array[0]
        );
        let end = format_args!(
            "flux(normalized)={:?}, θ={}",
            self.flux_array.get(last_index).unwrap() / flux_last,
            self.theta_array[last_index]
        );
        f.debug_struct("StationaryCurveSegment")
            .field("flux kind", &self.flux_kind())
            .field("start", &start)
            .field("end", &end)
            .field("length", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use dexter_machine::*;

    #[test]
    fn parametrized_segment_t() {
        let qfactor = UnityQfactor::new(LastClosedFluxSurface::Toroidal(0.1));
        let current = LarCurrent::new();
        let bfield = LarBfield::new();
        let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();

        let flux = MagneticFluxArray::Toroidal(Array1::linspace(1e-5, 0.09, 5));
        let theta = Array1::linspace(0.0, 3.1415, 5);
        let mut segment = ParametrizedSegment::new(flux, theta);

        assert_eq!(segment.len(), 5);

        let expected_t = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(segment.t.abs_diff_eq(&expected_t, 1e-14));

        segment.calculate_theta_dot(machine, -0.0123, 1e-5).unwrap();

        assert!(
            segment
                .theta_dot_array()
                .unwrap()
                .iter()
                .all(|theta_dot| theta_dot.is_finite())
        )
    }
}
