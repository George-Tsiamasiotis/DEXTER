//! Defines [`ParametrizedSegment`], a container type to store and parametrize the Stationary
//! Curve's distinct segments.
//!
//! [`StationaryCurveSegment`] is built on top of [`ParametrizedSegment`].

use std::f64::consts::TAU;

use crate::{CriticalPoint, GCHamiltonian};
use dexter_machine::{Accelerator, Accelerator2d, Machine, MagneticFluxArray, MagneticFluxKind};
use ndarray::{Array1, ArrayRef1};
use parabola::{Line, Point};
use roots::{Convergency, SearchError, find_root_brent, find_root_secant};
use rsl_interpolation::{Interpolation, LinearInterpolator};

/// Convergency criteria for better fitting the segments produced by the contour generator on the
/// true stationary curve.
struct ContourConvergency;

impl ContourConvergency {
    /// Maximum allowed value for `abs(𝜕B/𝜕𝜃)`.
    const EPS: f64 = 1e-16;
    /// Maximum allowed iterations.
    const MAX_ITER: usize = 10;
    /// Maximum allowed difference between the last elements `theta[i]` and `theta[i-1]`.
    ///
    /// Since the parameter `t` correspond to the array indices, a value of `1e-6` means that the
    /// the series converges when `ψ(t)` and `θ(t)` are approximated within `.0001%` of the initial
    /// spacing.
    const MAX_PARAMETER_DIFF: f64 = 1e-8;
}

#[expect(clippy::renamed_function_params, reason = "clarity")]
impl Convergency<f64> for ContourConvergency {
    fn is_root_found(&mut self, db_dtheta: f64) -> bool {
        db_dtheta.abs() < Self::EPS.abs()
    }
    fn is_converged(&mut self, theta_lo: f64, theta_hi: f64) -> bool {
        (theta_lo - theta_hi).abs() < Self::MAX_PARAMETER_DIFF
    }
    fn is_iteration_limit_reached(&mut self, iter: usize) -> bool {
        iter >= Self::MAX_ITER
    }
}

/// Helper struct to keep track of the errors encountered during the
/// [`ParametrizedSegment::improve`] method.
#[derive(Default, Clone, Copy)]
struct ImproveStats {
    /// From [`SearchError::NoConvergency`].
    no_convergency: usize,
    /// From [`SearchError::NoBracketing`].
    no_bracketing: usize,
    /// From [`SearchError::ZeroDerivative`].
    zero_derivative: usize,
    /// Average number of steps over all analyses.
    average_steps: f64,
}

// ===============================================================================================

/// Convergency criteria for finding the root of `dot(θ)(t)`.
struct FixedPointConvergency;

impl FixedPointConvergency {
    /// Maximum allowed value for `abs(dot(θ))`.
    const EPS: f64 = 1e-16;
    /// Maximum allowed iterations.
    const MAX_ITER: usize = 50;
    /// Maximum allowed difference between the 2 last `t` values.
    const MAX_PARAMETER_DIFF: f64 = 1e-7;
}

#[expect(clippy::renamed_function_params, reason = "clarity")]
impl Convergency<f64> for FixedPointConvergency {
    fn is_root_found(&mut self, theta_dot: f64) -> bool {
        theta_dot.abs() < Self::EPS.abs()
    }
    fn is_converged(&mut self, tlo: f64, thi: f64) -> bool {
        (tlo - thi).abs() < Self::MAX_PARAMETER_DIFF
    }
    fn is_iteration_limit_reached(&mut self, iter: usize) -> bool {
        iter >= Self::MAX_ITER
    }
}

// ===============================================================================================

/// A 2D parametrized curve on the (ψ, θ) poloidal plane.
///
/// The parameter `t` is chosen to be the index of the corresponding `(ψ,θ)` point.
#[derive(Clone)]
pub(crate) struct ParametrizedSegment {
    /// The segment's parameter.
    t: Array1<f64>,
    /// The segment's magnetic flux values.
    flux_array: MagneticFluxArray,
    /// The segment's `θ` values.
    theta_array: Array1<f64>,
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
}

// ===============================================================================================

/// A single, continuous segment of the [`StationaryCurve`][crate::StationaryCurve].
#[derive(Debug, Clone)]
pub struct StationaryCurveSegment {
    /// The (t, ψ, θ) parametrized segment.
    segment: ParametrizedSegment,
    /// Statistics of the segment's fitting on the true stationary curve.
    improve_stats: ImproveStats,
    /// The residual sum of squares of the `𝜕B/𝜕𝜃` values. The smaller the value, the better the fitting.
    rss: f64,
}

impl StationaryCurveSegment {
    /// Creates a new `StationaryCurveSegment` from the `flux` and `theta` arrays.
    ///
    /// # Panics
    ///
    /// This method panics if the input arrays do not have the same length or have 0 elements.
    #[must_use]
    pub(crate) fn new(flux_array: MagneticFluxArray, theta_array: Array1<f64>) -> Self {
        let segment = ParametrizedSegment::new(flux_array, theta_array);
        let improve_stats = ImproveStats::default();
        let rss = f64::NAN;
        Self {
            segment,
            improve_stats,
            rss,
        }
    }

    /// Returns the length of the segment.
    #[must_use]
    #[expect(clippy::len_without_is_empty, reason = "cannot create empty segments")]
    pub fn len(&self) -> usize {
        self.segment.len()
    }

    /// Returns the kind of the contained `flux_array`.
    #[must_use]
    pub fn flux_kind(&self) -> MagneticFluxKind {
        self.segment.flux_kind()
    }

    /// Returns a reference to the contained `flux_array`.
    #[must_use]
    pub fn flux_array(&self) -> &MagneticFluxArray {
        &self.segment.flux_array
    }

    /// Returns a reference to the contained `theta_array`.
    #[must_use]
    pub fn theta_array(&self) -> &Array1<f64> {
        &self.segment.theta_array
    }
}

/// Calculations.
impl StationaryCurveSegment {
    /// Calculates the segment's critical points, that is, the points for which `dot(θ)(t) = 0`.
    ///
    /// Note that no intermediate steps are stored.
    #[must_use]
    pub fn calculate_critical_points(
        &self,
        machine: Machine,
        pzeta: f64,
        mu: f64,
    ) -> Vec<CriticalPoint> {
        let theta_dot_array = self.calculate_theta_dot(machine, pzeta, mu);
        let root_indices = Self::find_root_indices(&theta_dot_array);
        let exact_roots = self.find_exact_roots(machine, root_indices, pzeta, mu);
        self._calculate_critical_points(machine, exact_roots, pzeta, mu)
    }

    /// Calculates the `dot(θ)` array from the contained `flux` and `theta` arrays and the constants
    /// of motion `Pζ` and `μ`.
    #[must_use]
    fn calculate_theta_dot(&self, machine: Machine, pzeta: f64, mu: f64) -> Array1<f64> {
        let Ok(theta_dot_array) = GCHamiltonian::eval_theta_dot_array(
            machine,
            &self.segment.theta_array,
            &self.segment.flux_array,
            pzeta,
            mu,
        ) else {
            unreachable!("By construction theta and flux are in-bounds and flux is good")
        };
        theta_dot_array
    }

    /// Brings each of the segment's points (as calculated by the contour generator) closer to
    /// the true stationary curve.
    ///
    /// This is achieved by finding the root of `𝜕B/𝜕𝜃 = 0` close to the initial point using an
    /// iterative method.
    ///
    /// The contour point and its two adjacent points form a triangle. It is highly likely that the
    /// true stationary curve passes through this triangle. We then restrict `ψ` and `θ` to move on
    /// the line defined by the initial point and the midpoint of the two adjacent points, thus
    /// reducing the problem to the 1D root finding problem `𝜕B/𝜕𝜃(ψ(θ), θ) = 0`.
    ///
    /// We use the `secant` method with initial points the point itself and the midpoint of its two
    /// adjacent contour points, which lie very close to the true curve. The secant method permits
    /// stepping outside this interval in case the true curve passes outside the triangle.
    pub(crate) fn improve(&mut self, machine: Machine) {
        let len = self.len();
        // skip the first and last points as we cannot define the triangle. These points always lie
        // on the axis or on the wall.
        (0..len).skip(1).take(len - 2).for_each(|i| {
            let x = self.segment.theta_array[i];
            let y = self.segment.flux_array[i];
            let xlo = self.segment.theta_array[i - 1];
            let xhi = self.segment.theta_array[i + 1];
            let ylo = self.segment.flux_array[i - 1];
            let yhi = self.segment.flux_array[i + 1];
            let xmid = xlo.midpoint(xhi);
            let ymid = ylo.midpoint(yhi);

            // Define `ψ(θ)` line
            if xmid == 0.0 {
                return;
            };
            let slope = (ymid - y) / (xmid - x);
            let point = Point { x: xmid, y: ymid };
            let line_equation = Line::from_slope_point(slope, point);
            let line = |theta| line_equation.eval(theta);

            // Define function over line
            let mut steps = 0.0;
            let acc = &mut Accelerator2d::new();
            let fun = |theta: f64| {
                let theta = theta.rem_euclid(TAU);
                let flux = self.flux_kind().to_magnetic_flux(line(theta));
                match machine.bfield().eval_deriv_theta(flux, theta, acc) {
                    Ok(value) => {
                        steps += 1.0;
                        value
                    }
                    Err(_) => f64::MAX, // Stops the root finder. Should be unreachable
                }
            };

            let mut convergency = ContourConvergency;
            let root = match find_root_secant(x, xmid, fun, &mut convergency) {
                Ok(root) => {
                    self.improve_stats.average_steps += steps / self.len() as f64;
                    root
                }
                Err(err) => {
                    match err {
                        SearchError::NoConvergency => self.improve_stats.no_convergency += 1,
                        SearchError::NoBracketing => self.improve_stats.no_bracketing += 1,
                        SearchError::ZeroDerivative => self.improve_stats.zero_derivative += 1,
                    }
                    return;
                }
            };

            let root_flux = line(root);
            self.segment.theta_array[i] = root.rem_euclid(TAU);
            self.segment.flux_array.set_value(i, root_flux);
        });
    }

    /// Calculates the residual sum of squares of `𝜕B/𝜕𝜃` over the stored `theta` and `flux` arrays.
    pub(crate) fn calculate_rss(&mut self, machine: Machine) {
        let acc = &mut Accelerator2d::new();
        self.rss = self
            .segment
            .flux_array
            .values()
            .iter()
            .zip(&self.segment.theta_array)
            .map(|(flux, theta)| {
                let flux = self.flux_kind().to_magnetic_flux(*flux);
                match machine.bfield().eval_deriv_theta(flux, *theta, acc) {
                    Ok(db_dtheta) => db_dtheta,
                    Err(_) => unreachable!("contours are always in-bounds"),
                }
            })
            .collect::<Array1<f64>>()
            .powi(2)
            .sum()
    }

    /// Calculates the floored indices of the `dot(θ)(t)` roots.
    ///
    /// The parameter `t` was defined to correspond to the index of each `(ψ, θ)` point. This means
    /// that if `i` is a floored index, then a root exists in the interval `i <= t <= i+1`.
    #[expect(clippy::pattern_type_mismatch, reason = "windows")]
    fn find_root_indices(theta_dot_array: &ArrayRef1<f64>) -> Box<[usize]> {
        let slice = theta_dot_array
            .as_slice()
            .expect("non-empty by construction");
        slice
            .windows(2)
            .enumerate()
            .filter(|(_, window)| window[0].signum() != window[1].signum())
            .map(|(i, _)| i)
            .collect()
    }

    /// Finds the true roots of `dot(θ)(t)`.
    ///
    /// For every `t[i]`, `t[i+1]` interval where `dot(θ)` changes sign, we construct the 1D function:
    /// `f(t) = f(ψ(t), θ(t)) = dot(θ)(t)` and find its root.
    pub(crate) fn find_exact_roots(
        &self,
        machine: Machine,
        root_indices: Box<[usize]>,
        pzeta: f64,
        mu: f64,
    ) -> Box<[f64]> {
        let mut exact_roots = Vec::new();
        root_indices.into_iter().for_each(|index| {
            // Index into the arrays at the edges of the `[i, i+1]` interval. By construction of
            // the parameter `t`, it always holds that `tlo<thi` and therefore the interpolation
            // is always defined.
            let tlo = self.segment.t[index]; // == index
            let thi = self.segment.t[index + 1]; // == index + 1
            let theta_lo = self.segment.theta_array[index];
            let theta_hi = self.segment.theta_array[index + 1];
            let flux_lo = self.segment.flux_array[index];
            let flux_hi = self.segment.flux_array[index + 1];
            let ta = [tlo, thi];

            let ya_theta = [theta_lo, theta_hi];
            let ya_flux = [flux_lo, flux_hi];

            #[expect(clippy::float_cmp, reason = "need exact equality")]
            if tlo == thi {
                // ZeroDerivative error
                return;
            }

            let acc = &mut Accelerator2d::new();

            let fun = |t: f64| {
                if !(tlo..=thi).contains(&t) {
                    return f64::MAX; // stop the solver. Landing here should be extremely unlikely
                }

                // `t -> θ` interpolation
                let theta = LinearInterpolator
                    .eval(&ta, &ya_theta, t, acc.xacc())
                    .expect("arrays and values ok");

                // `t -> ψ/ψp` interpolation
                let flux = LinearInterpolator
                    .eval(&ta, &ya_flux, t, acc.xacc())
                    .expect("arrays and values ok");

                let flux = self.flux_kind().to_magnetic_flux(flux);
                GCHamiltonian::eval_theta_dot(machine, theta, flux, pzeta, mu).unwrap_or(f64::NAN)
            };

            let mut convergency = FixedPointConvergency;
            let Ok(t_root) = find_root_brent(tlo, thi, fun, &mut convergency) else {
                return;
            };

            exact_roots.push(t_root);
        });
        exact_roots.into_boxed_slice()
    }

    /// Calculates the critical points by interpolating each exact root to the corresponding `θ` and
    /// `ψ/ψp` points.
    pub(crate) fn _calculate_critical_points(
        &self,
        machine: Machine,
        exact_roots: Box<[f64]>,
        pzeta: f64,
        mu: f64,
    ) -> Vec<CriticalPoint> {
        exact_roots
            .into_iter()
            .filter_map(|root| {
                // Index into the arrays at the edges of the `[i, i+1]` interval. By construction of
                // the parameter `t`, it always holds that `tlo<thi` and therefore the interpolation
                // is always defined.
                let floor_index = root.floor() as usize;
                let ceil_index = root.ceil() as usize;
                let tlo = self.segment.t[floor_index];
                let thi = self.segment.t[ceil_index];
                let theta_lo = self.segment.theta_array[floor_index];
                let theta_hi = self.segment.theta_array[ceil_index];
                let flux_lo = self.segment.flux_array.values()[floor_index];
                let flux_hi = self.segment.flux_array.values()[ceil_index];
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

                let flux = self.segment.flux_array.kind().to_magnetic_flux(flux);
                CriticalPoint::new(machine, theta, flux, pzeta, mu).ok()
            })
            .collect()
    }
}

impl std::fmt::Debug for ImproveStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ImproveStats")
            .field("no_convergency", &self.no_convergency)
            .field("no_bracketing", &self.no_bracketing)
            .field("zero_derivative", &self.zero_derivative)
            .field("average_steps", &format!("{:.5}", self.average_steps))
            .finish()
    }
}

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
        let found_indices = StationaryCurveSegment::find_root_indices(&theta_dot_array);
        assert_eq!(&found_indices, &expected_indices);
    }
}
