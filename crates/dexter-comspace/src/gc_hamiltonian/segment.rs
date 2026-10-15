//! Defines [`ParametrizedSegment`], a container type to store and parametrize the Stationary
//! Curve's distinct segments.

use dexter_machine::{Accelerator, Machine, MagneticFluxArray, MagneticFluxKind};
use ndarray::{Array1, ArrayView1};
use rsl_interpolation::{BuildInterpolator, Interpolation, LinearInterpolator};

use crate::{CriticalPoint, GCHamiltonian};

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

    /// Calculates the segment's, critical points, that is, the points for which `dot(θ)(t) = 0`.
    ///
    /// Note that no intermediate steps are stored other than `dot(θ)(t)` as evaluated over all
    /// `(θ, ψ/ψp)` pairs.
    pub fn calculate_critical_points(
        &mut self,
        machine: Machine,
        pzeta: f64,
        mu: f64,
    ) -> Vec<CriticalPoint> {
        self.calculate_theta_dot(machine, pzeta, mu);
        let root_indices = self.find_root_indices();
        let exact_roots = self.find_exact_roots(root_indices);
        self._calculate_critical_points(machine, pzeta, mu, exact_roots)
    }
}

/// Calculations.
impl ParametrizedSegment {
    /// Calculates the `dot(θ)` array from the contained `flux` and `theta` arrays and the constants
    /// of motion `Pζ` and `μ`.
    pub fn calculate_theta_dot(&mut self, machine: Machine, pzeta: f64, mu: f64) {
        let Ok(theta_dot_array) = GCHamiltonian::eval_theta_dot_array(
            machine,
            &self.theta_array,
            &self.flux_array,
            pzeta,
            mu,
        ) else {
            unreachable!("By construction theta and flux are in-bounds and flux is good")
        };
        self.theta_dot = Some(theta_dot_array);
    }

    /// Calculates the floored indices of the `dot(θ)(t)` roots.
    ///
    /// The parameter `t` was defined to correspond to the index of each `(ψ, θ)` point. This means
    /// that if `i` is a floored index, then a root exists in the interval `i <= t <= i+1`.
    #[expect(clippy::pattern_type_mismatch, reason = "windows")]
    fn find_root_indices(&self) -> Box<[usize]> {
        let slice = self
            .theta_dot
            .as_ref()
            .expect("method called after theta dot calculation")
            .as_slice()
            .expect("non-empty by construction");
        slice
            .windows(2)
            .enumerate()
            .filter(|(_, window)| window[0].signum() != window[1].signum())
            .map(|(i, _)| i)
            .collect()
    }

    /// Finds the true roots of `dot(θ)(t)` by interpolating `t` on the `[i, i+1]` intervals, where
    /// `i` the already calculated floored indices.
    ///
    /// The interpolator is a linear interpolator defined as `dot(θ) -> t` within the `[i, i+1]`
    /// interval. and calculates `t(dot(θ) = 0)`.
    pub(crate) fn find_exact_roots(&self, root_indices: Box<[usize]>) -> Box<[f64]> {
        use std::cmp::Ordering::*;
        let theta_dot = self
            .theta_dot
            .as_ref()
            .expect("method called after theta dot calculation");

        root_indices
            .into_iter()
            .filter_map(|index| -> Option<f64> {
                // Define the `[x0, x1] -> [y0, y1]` intervals.
                let xlo = theta_dot[index];
                let xhi = theta_dot[index + 1];
                let ylo = self.t[index];
                let yhi = self.t[index + 1];

                // In the case where `x0>x1`, reverse the intervals. This does not affect the final
                // interpolated value.
                //
                // In the very edge case where `xlo=xhi`, the interpolator might either panic or
                // return garbage, so we simply discard the root.
                let (xa, ya) = match xlo.total_cmp(&xhi) {
                    Less => ([xlo, xhi], [ylo, yhi]),
                    Greater => ([xhi, xlo], [yhi, ylo]),
                    Equal => return None,
                };

                let interp = LinearInterpolator::build(&xa, &ya).expect("arrays are ok");
                interp.eval(&xa, &ya, 0.0, &mut Accelerator::new()).ok()
            })
            .collect()
    }

    /// Calculates the critical points by interpolating each exact root to the corresponding `θ` and
    /// `ψ/ψp` points.
    pub(crate) fn _calculate_critical_points(
        &self,
        machine: Machine,
        pzeta: f64,
        mu: f64,
        exact_roots: Box<[f64]>,
    ) -> Vec<CriticalPoint> {
        exact_roots
            .into_iter()
            .filter_map(|root| {
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
                    .eval(&xa, &ya_theta, root, acc)
                    .expect("arrays and values ok");

                // `t -> ψ/ψp` interpolation
                let ya_flux = [flux_lo, flux_hi];
                let flux = LinearInterpolator
                    .eval(&xa, &ya_flux, root, acc)
                    .expect("arrays and values ok");

                let flux = self.flux_array.kind().to_magnetic_flux(flux);
                CriticalPoint::new(machine, theta, flux, pzeta, mu).ok()
            })
            .collect()
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
    use ndarray::array;

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

        segment.calculate_theta_dot(machine, -0.0123, 1e-5);

        assert!(
            segment
                .theta_dot_array()
                .unwrap()
                .iter()
                .all(|theta_dot| theta_dot.is_finite())
        )
    }

    #[test]
    fn find_roots() {
        let flux = Array1::linspace(0.0, 1.0, 8);
        let theta = Array1::linspace(1.0, 2.0, 8);
        let theta_dot = array![2.0, 1.0, -1.0, -0.5, 0.3, 1.0, 2.0, -3.0];
        let mut segment = ParametrizedSegment::new(MagneticFluxArray::Toroidal(flux), theta);
        segment.theta_dot = Some(theta_dot);

        let expected_indices: Box<[usize]> = Box::new([1, 3, 6]);
        let found_indices = segment.find_root_indices();
        assert_eq!(&found_indices, &expected_indices);

        let _expected_roots = array![0.25, 0.0, -0.5]; // roots of `dot(θ)(t)`
        let expected_roots_t = array![1.0 + (1.0 / 2.0), 3.0 + (0.5 / 0.8), 6.0 + (2.0 / 5.0)];

        let found_roots = Array1::from(segment.find_exact_roots(found_indices));

        assert!(found_roots.relative_eq(&expected_roots_t, 1e-12, 1e-12))
    }
}
