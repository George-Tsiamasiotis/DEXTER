//! Defines [`ParametrizedSegment`], a container type to store and parametrize 2D curves on the
//! (θ, ψ) poloidal plane.

use dexter_machine::{Accelerator, MagneticFlux, MagneticFluxArray, MagneticFluxKind};
use ndarray::{Array1, ArrayRef1};
use rsl_interpolation::{Interpolation, LinearInterpolator};

// ===============================================================================================

/// A 2D parametrized curve on the (ψ, θ) poloidal plane.
///
/// The parameter `t` is chosen to be the index of the corresponding `(ψ,θ)` point.
#[derive(Clone)]
pub(crate) struct ParametrizedSegment {
    /// The segment's parameter.
    pub(crate) t: Array1<f64>,
    /// The segment's magnetic flux values.
    pub(crate) flux_array: MagneticFluxArray,
    /// The segment's `θ` values.
    pub(crate) theta_array: Array1<f64>,
}

impl ParametrizedSegment {
    /// Creates a new `ParametrizedSegment` from the `flux` and `theta` arrays and initializes the
    /// parameter `t`.
    ///
    /// # Panics
    ///
    /// This method panics if the input arrays do not have the same length or have 0 elements.
    #[must_use]
    pub(crate) fn new(flux_array: MagneticFluxArray, theta_array: Array1<f64>) -> Self {
        assert_eq!(
            flux_array.len(),
            theta_array.len(),
            "Input arrays must have the same length"
        );
        assert_ne!(theta_array.len(), 0, "Input arrays cannot be empty");
        let t = Array1::linspace(0.0, (flux_array.len() - 1) as f64, flux_array.len());
        Self {
            t,
            flux_array,
            theta_array,
        }
    }

    /// Returns the number of points that comprise the segment.
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.t.len()
    }

    /// Returns the [`MagneticFluxArray`] of the stored [`MagneticFluxArray`].
    #[must_use]
    pub(crate) fn flux_kind(&self) -> MagneticFluxKind {
        self.flux_array.kind()
    }

    /// Calculates the floored indices of the `t` array where the array changes sign.
    ///
    /// The parameter `t` was defined to correspond to the index of each `(ψ, θ)` point. This means
    /// that if `i` is a floored index, then a root exists in the interval `i <= t <= i+1`.
    #[expect(clippy::pattern_type_mismatch, reason = "windows")]
    pub(crate) fn find_root_indices(array: &ArrayRef1<f64>) -> Box<[usize]> {
        let slice = array.as_slice().expect("non-empty by construction");
        slice
            .windows(2)
            .enumerate()
            .filter(|(_, window)| window[0].signum() != window[1].signum())
            .map(|(i, _)| i)
            .collect()
    }

    /// Performs the `t -> (θ, ψ)` interpolation in every discorvered exact root.
    pub(crate) fn interpolate_exact_roots(&self, exact_roots: &[f64]) -> Vec<(f64, MagneticFlux)> {
        exact_roots
            .iter()
            .map(|root| {
                // Index into the arrays at the edges of the `[i, i+1]` interval. By construction of
                // the parameter `t`, it always holds that `tlo<thi` and therefore the interpolation
                // is always defined.
                let floor_index = root.floor() as usize;
                let ceil_index = root.ceil() as usize;
                let tlo = self.t[floor_index];
                let thi = self.t[ceil_index];
                let theta_lo = self.theta_array[floor_index];
                let theta_hi = self.theta_array[ceil_index];
                let flux_lo = self.flux_array.values()[floor_index];
                let flux_hi = self.flux_array.values()[ceil_index];
                let xa = [tlo, thi];

                let acc = &mut Accelerator::new(); // doesn't really do anything

                // `t -> θ` interpolation
                let ya_theta = [theta_lo, theta_hi];
                let theta = LinearInterpolator
                    .eval(&xa, &ya_theta, *root, acc)
                    .expect("arrays and values ok");

                // `t -> ψ` interpolation
                let ya_flux = [flux_lo, flux_hi];
                let flux = LinearInterpolator
                    .eval(&xa, &ya_flux, *root, acc)
                    .expect("arrays and values ok");

                let flux = self.flux_array.kind().to_magnetic_flux(flux);
                (theta, flux)
            })
            .collect()
    }
}

// ===============================================================================================

impl std::fmt::Debug for ParametrizedSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let last_index = self.len() - 1;
        let start = format_args!(
            "flux={}: {:.8}, θ={:.8}",
            self.flux_kind().symbol(),
            self.flux_array[0],
            self.theta_array[0]
        );
        let end = format_args!(
            "flux={}: {:.8}, θ={:.8}",
            self.flux_kind().symbol(),
            self.flux_array[last_index],
            self.theta_array[last_index]
        );
        f.debug_struct("StationaryCurveSegment")
            .field("start", &start)
            .field("end", &end)
            .field("length", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ndarray::array;

    #[test]
    fn parametrized_segment_t() {
        let flux = MagneticFluxArray::Toroidal(Array1::linspace(1e-5, 0.09, 5));
        let theta = Array1::linspace(0.0, 3.1415, 5);
        let segment = ParametrizedSegment::new(flux, theta);

        assert_eq!(segment.len(), 5);

        let expected_t = Array1::from_vec(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(segment.t.abs_diff_eq(&expected_t, 1e-14));
    }

    #[test]
    fn find_root_indices() {
        let theta_dot_array = array![2.0, 1.0, -1.0, -0.5, 0.3, 1.0, 2.0, -3.0];

        let expected_indices: Box<[usize]> = Box::new([1, 3, 6]);
        let found_indices = ParametrizedSegment::find_root_indices(&theta_dot_array);
        assert_eq!(&found_indices, &expected_indices);
    }
}
