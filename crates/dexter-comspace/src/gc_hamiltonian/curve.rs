//! Defines the [`StationaryCurve`] and [`StationaryCurveSegment`] types.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;
use std::ops::Range;
use std::time::{Duration, Instant};

use contour::ContourBuilder;
use dexter_machine::{Accelerator2d, Machine, MagneticFluxArray, MagneticFluxKind};
use geo_types::{LineString, MultiLineString};
use ndarray::{Array1, Array2, ArrayRef1};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rsl_interpolation::{Interpolation, LinearInterpolator};

use crate::gc_hamiltonian::ParametrizedSegment;
use crate::{CriticalPoint, GCHamiltonian};
use roots::{Convergency, find_root_brent, find_root_secant};

/// The density of the (θ, ψ/ψp) grid on which to search the stationary curve contours.
const SC_CONTOUR_DENSITY: (usize, usize) = (300, 300);

// ===============================================================================================

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

// ===============================================================================================

/// Convergency criteria for finding the root of `dot(θ)(t)`.
struct FixedPointConvergency;

impl FixedPointConvergency {
    /// Maximum allowed value for `abs(dot(θ))`.
    const EPS: f64 = 1e-16;
    /// Maximum allowed iterations.
    const MAX_ITER: usize = 50;
    /// Maximum allowed difference between the 2 last `t` values.
    const MAX_PARAMETER_DIFF: f64 = 1e-9;
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

/// The stationary curve of the Hamiltonian.
///
/// The stationary curve is defined through the equation `𝜕H/𝜕𝜃 = 0`. In the absence of an electric
/// field, this is equivalent to `𝜕B/𝜕𝜃 = 0`.
///
/// The stationary curve is described by one or more [`segments`][StationaryCurveSegment]
/// of the form `𝜕H/𝜕𝜃(ψ,θ)=0`.
#[derive(Clone)]
pub struct StationaryCurve {
    /// The separated segments comprising the `StationaryCurve`.
    segments: Box<[StationaryCurveSegment]>,
    /// The duration of the [`Self::calculate_critical_points`] routine.
    duration: Option<Duration>,
}

impl StationaryCurve {
    /// Creates a new `StationaryCurve`.
    ///
    /// # Process
    ///
    /// 1. Build the (θ, ψ/ψp) *orthogonal* grid.
    /// 2. Set up the [`ContourBuilder`].
    /// 3. Extract [`ContourBuilder::lines`] into a [`Box<StationaryCurveSegment>`] and store it.
    ///
    /// # Panics
    ///
    /// This method panics if [`ContourBuilder::lines`] returns an error. This is unusual since at
    /// least one segment of `𝜕B/𝜕𝜃 = 0` should always exist and therefore should be fatal.
    #[must_use]
    pub fn new(machine: Machine) -> Self {
        Self::_build(machine)
    }

    /// Returns a reference to the stored [`StationaryCurveSegment`]s.
    #[must_use]
    pub fn segments(&self) -> &[StationaryCurveSegment] {
        &self.segments
    }

    /// Returns the number of stored [`StationaryCurveSegment`]s.
    #[must_use]
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Returns `true` if the curve contains no [`StationaryCurveSegment`]s.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.segments.len() == 0
    }

    /// Returns the duration of the [`Self::calculate_critical_points`] routine.
    #[must_use]
    pub fn duration(&self) -> Option<Duration> {
        self.duration
    }

    /// Calculates the critical points of the `(E, Pζ; μ)` slice over all values of `pzeta_array` and
    /// stores them.
    pub fn calculate_critical_points(
        &mut self,
        machine: Machine,
        pzeta_array: &ArrayRef1<f64>,
        mu: f64,
    ) -> Box<[CriticalPoint]> {
        let start = Instant::now();

        let slice = pzeta_array.as_slice().unwrap_or(&[]); // do nothing
        let points: Vec<CriticalPoint> = slice
            .into_par_iter()
            .map(|pzeta| {
                self.segments
                    .par_iter()
                    .map(|segment| segment.calculate_critical_points(machine, *pzeta, mu))
                    .flatten()
                    .collect::<Vec<CriticalPoint>>()
            })
            .flatten()
            .collect();

        self.duration = Some(start.elapsed());

        points.into_boxed_slice()
    }

    /// Builds the `StationaryCurve` of `machine`.
    pub(crate) fn _build(machine: Machine) -> Self {
        let (theta_array, flux_array, grid) = build_grid(machine);
        let flux_arrayv = flux_array.values();
        let flux_cutoff = flux_arrayv.last().copied().expect("non-empty");

        let dx = theta_array[1] - theta_array[0];
        let dy = flux_arrayv[1] - flux_arrayv[0];

        // NOTE: `ContourBuilder` defines the position of each pixel as its center, therefore
        // we must shift it back by half a step in both dimensions.
        let x0 = theta_array[0] - dx / 2.0;
        let y0 = flux_arrayv[0] - dy / 2.0;

        let cb = ContourBuilder::new(theta_array.len(), flux_array.len(), true)
            .x_origin(x0)
            .y_origin(y0)
            .x_step(dx)
            .y_step(dy);

        let Some(values) = grid.as_slice() else {
            unreachable!("contiguous and standard layout by construction (ndarray default)")
        };

        let isolines = match cb.lines(values, &[0.0]) {
            Ok(isolines) => isolines,
            Err(err) => panic!("ContourBuilder error: {err}"),
        };

        // Post-processing of contour lines into `ParametrizedSegments`
        let strings = isolines_to_strings(&isolines);
        let mut separated_strings = strings_into_separated_strings(strings, flux_cutoff);
        clamp_separated_strings(&mut separated_strings);
        let mut segments = separated_strings_into_segments(separated_strings, machine.good_flux());

        segments
            .iter_mut()
            .for_each(|segment| segment.improve(machine));
        segments
            .iter_mut()
            .for_each(|segment| segment.calculate_rss(machine));

        Self {
            segments,
            duration: None,
        }
    }
}

/// Builds the (θ, ψ/ψp) grid on which to search for stationary curve contours.
fn build_grid(machine: Machine) -> (Array1<f64>, MagneticFluxArray, Array2<f64>) {
    let flux_kind = machine.good_flux();
    let flux_last = machine.last(flux_kind);

    // Avoid edge phenomena
    let theta_first = -0.1;
    let theta_last = TAU + 0.1;
    let flux_first = 1e-9;
    let flux_last = (1.0 - 1e-9) * flux_last.value();

    // Arrays must be linspaces
    let theta_array = Array1::linspace(theta_first, theta_last, SC_CONTOUR_DENSITY.0);
    let flux_values = Array1::linspace(flux_first, flux_last, SC_CONTOUR_DENSITY.1);
    let flux_array = flux_kind.to_magnetic_flux_array(flux_values);

    let mut grid = Array2::<f64>::uninit((flux_array.len(), theta_array.len()));

    assert_eq!(grid.ncols(), theta_array.len(), "sanity check");
    assert_eq!(grid.nrows(), flux_array.len(), "sanity check");

    let acc = &mut Accelerator2d::new();

    for m in 0..grid.nrows() {
        let Some(flux) = flux_array.get(m) else {
            unreachable!("'m' is in-bounds by construction")
        };
        for n in 0..grid.ncols() {
            let theta = theta_array[n].rem_euclid(TAU);
            if let Ok(db_dtheta) = machine.bfield().eval_deriv_theta(flux, theta, acc) {
                grid[[m, n]] = MaybeUninit::new(db_dtheta)
            } else {
                unreachable!("arrays are in-bounds and flux is Good")
            }
        }
    }

    // SAFETY: The loop passes from all elements and initializes them
    let grid = unsafe { grid.assume_init() };

    (theta_array, flux_array, grid)
}

/// Converts a slice of `Lines` into a `Box` of `LineStrings` without cloning.
fn isolines_to_strings(isolines: &[contour::Line]) -> Box<[&LineString]> {
    isolines
        .iter()
        .map(contour::Line::geometry)
        .collect::<Vec<&MultiLineString>>()
        .into_iter()
        .flatten()
        .collect()
}

/// Consumes a `Box` of `LineStrings`, passes every line string through [`separate`] and returns the
/// new separated segments.
fn strings_into_separated_strings(
    strings: Box<[&LineString]>,
    flux_cutoff: f64,
) -> Box<[LineString]> {
    strings
        .into_iter()
        .map(|string| separate(string, flux_cutoff))
        .collect::<Vec<Box<[LineString]>>>()
        .into_iter()
        .flatten()
        .collect()
}

/// Iterates through each `separated_string`, clamping any slightly negative `theta` that may appear
/// to 0.0.
///
/// Slightly negative `thetas` can appear when a stationary curve segment lies exactly on the `θ=0`
/// line, which is often the case.
fn clamp_separated_strings(separated_strings: &mut [LineString]) {
    let clamp_range = -1e-4 * TAU..0.0;

    separated_strings.iter_mut().for_each(|line_string| {
        line_string.coords_mut().for_each(|coord| {
            if clamp_range.contains(&coord.x) {
                coord.x = 0.0
            }
        });
    });
}

/// Consumes a `Box` of separated `LineStrings` and maps each line string into a
/// [`ParametrizedSegment`].
fn separated_strings_into_segments(
    separated_strings: Box<[LineString]>,
    flux_kind: MagneticFluxKind,
) -> Box<[StationaryCurveSegment]> {
    separated_strings
        .into_iter()
        .map(|separated_string| {
            let len = separated_string.0.len();
            let mut theta = Vec::<f64>::with_capacity(len);
            let mut flux = Vec::<f64>::with_capacity(len);

            separated_string.coords().for_each(|coord| {
                theta.push(coord.x);
                flux.push(coord.y);
            });

            let theta_array = Array1::from_vec(theta);
            let flux_arrayv = Array1::from_vec(flux);
            let flux_array = flux_kind.to_magnetic_flux_array(flux_arrayv);

            StationaryCurveSegment::new(flux_array, theta_array)
        })
        .collect()
}

/// Separates a [`LineString`] into multiple `LineStrings` by removing the parts that touch
/// the wall and lie outside the [0,2π] interval.
///
/// This is necessary as [`contour`] only yields closed surfaces. However, the "useless"
/// segments we want to discard always connect the true isolines by hugging the grid's
/// bounds (By construction of the contour lines through the marching squares algorithm,
/// these "connecting" segments cannot lie anywhere but on the grid's bounds.). We can
/// use this to iterate over every string and safely separate each segment.
///
/// Any new lines with a length less that 1 are discarded.
fn separate(string: &LineString, flux_last: f64) -> Box<[LineString]> {
    let mut res = Vec::<LineString>::new();

    let xbound: Range<f64> = (-1e-4 * TAU)..((1.0 - 1e-4) * TAU);
    let ybound: Range<f64> = (1e-8 * flux_last)..((1.0 - 1e-10) * flux_last);

    // Keeps track of the end index of every discovered valid segment
    let mut end = 0;
    while end < string.0.len() {
        let iterator = string
            .coords()
            .skip(end)
            .take_while(|coord| {
                end += 1;
                xbound.contains(&coord.x) && ybound.contains(&coord.y)
            })
            .copied();
        let new_line = LineString::from_iter(iterator);
        if new_line.0.len() > 1 {
            res.push(new_line);
        }
    }
    res.into_boxed_slice()
}

// ===============================================================================================

/// A single, continuous segment of the [`StationaryCurve`][crate::StationaryCurve].
#[derive(Debug, Clone)]
pub struct StationaryCurveSegment {
    /// The (t, ψ, θ) parametrized segment.
    segment: ParametrizedSegment,
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
        let rss = f64::NAN;
        Self { segment, rss }
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
        let root_indices = ParametrizedSegment::find_root_indices(&theta_dot_array);
        let exact_roots = self.find_exact_roots(machine, root_indices, pzeta, mu);
        self.critical_points_from_exact_roots(machine, &exact_roots, pzeta, mu)
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

            let Ok(t_root) = find_root_brent(tlo, thi, fun, &mut FixedPointConvergency) else {
                return;
            };

            exact_roots.push(t_root);
        });
        exact_roots.into_boxed_slice()
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
            let point = parabola::Point { x: xmid, y: ymid };
            let line_equation = parabola::Line::from_slope_point(slope, point);
            let line = |theta| line_equation.eval(theta);

            // Define function over line
            let acc = &mut Accelerator2d::new();
            let fun = |theta: f64| {
                let theta = theta.rem_euclid(TAU);
                let flux = self.flux_kind().to_magnetic_flux(line(theta));
                machine
                    .bfield()
                    .eval_deriv_theta(flux, theta, acc)
                    .unwrap_or(f64::NAN) // stop the solver
            };

            let mut convergency = ContourConvergency;
            let Ok(root) = find_root_secant(x, xmid, fun, &mut convergency) else {
                return;
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

    /// Calculates the critical points by interpolating each exact root to the corresponding `θ` and
    /// `ψ/ψp` points.
    pub(crate) fn critical_points_from_exact_roots(
        &self,
        machine: Machine,
        exact_roots: &[f64],
        pzeta: f64,
        mu: f64,
    ) -> Vec<CriticalPoint> {
        self.segment
            .interpolate_exact_roots(exact_roots)
            .into_iter()
            .filter_map(|(theta, flux)| CriticalPoint::new(machine, theta, flux, pzeta, mu).ok())
            .collect()
    }
}

// ===============================================================================================

impl std::fmt::Debug for StationaryCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let duration_str: String = self
            .duration
            .map_or_else(|| "None".to_owned(), |d| format!("{d:?}"));
        f.debug_struct("StationaryCurve")
            .field(
                "flux_kind",
                &self.segments.first().map(StationaryCurveSegment::flux_kind),
            )
            .field("number of segments", &self.segments.len())
            .field("segments", &self.segments)
            .field("duration", &duration_str)
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use geo_types::line_string;

    const LAST: f64 = 10.0;

    fn assert_separated_equals_expected(string: LineString, expected: LineString) {
        let separated_strings = separate(&string, LAST);
        assert_eq!(separated_strings.len(), 1);
        let separated = separated_strings[0].clone();
        assert_eq!(separated, expected);
    }

    #[test]
    fn separate_line_string_none1() {
        let string: LineString = line_string![
            (x: -1.0, y:0.),
            (x: -1.0, y:0.),
        ];
        let separateped_strings = separate(&string, LAST);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none2() {
        let string: LineString = line_string![
            (x: -1., y:2000000.0),
            (x: -1., y:3000000.0),
        ];
        let separateped_strings = separate(&string, LAST);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none3() {
        let string: LineString = line_string![
            (x: -1., y:3.0),
        ];
        let separateped_strings = separate(&string, LAST);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none4() {
        let string: LineString = line_string![
            (x: -1., y:3.0),
            (x: -1., y:0.),
        ];
        let separateped_strings = separate(&string, LAST);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_single1() {
        let string: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        let expected: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        assert_separated_equals_expected(string, expected);
    }

    #[test]
    fn separate_line_string_single2() {
        let string: LineString = line_string![
            (x: 0., y:0.),
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        let expected: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        assert_separated_equals_expected(string, expected);
    }

    #[test]
    fn separate_line_string_single3() {
        let string: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:0.),
            (x: 0., y:0.),
        ];
        let expected: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        assert_separated_equals_expected(string, expected);
    }

    #[test]
    fn separate_line_string_single4() {
        let string: LineString = line_string![
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:0.),
        ];
        let expected: LineString = line_string![
            (x: 0., y:1.0),
            (x: 0., y:2.0),
        ];
        assert_separated_equals_expected(string, expected);
    }

    #[test]
    fn separate_line_string_multiple1() {
        let string: LineString = line_string![
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:0.),
            (x: 0., y:3.0),
            (x: 0., y:4.0),
            (x: 0., y:5.0),
            (x: 0., y:0.),
        ];
        let separated_strings = separate(&string, LAST);
        assert_eq!(separated_strings.len(), 2);

        let separated1 = separated_strings[0].clone();
        assert_eq!(
            separated1,
            line_string![
                (x: 0., y:1.0),
                (x: 0., y:2.0),
            ]
        );

        let separated2 = separated_strings[1].clone();
        assert_eq!(
            separated2,
            line_string![
                (x: 0., y:3.0),
                (x: 0., y:4.0),
                (x: 0., y:5.0),
            ]
        );
    }

    #[test]
    fn separate_line_string_multiple2() {
        let string: LineString = line_string![
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:1000000.0),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
            (x: 0., y:4.0),
            (x: 0., y:1000000.0),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
        ];
        let separated_strings = separate(&string, LAST);
        assert_eq!(separated_strings.len(), 4);

        let separated1 = separated_strings[0].clone();
        assert_eq!(
            separated1,
            line_string![
                (x: 0., y:1.0),
                (x: 0., y:2.0),
            ]
        );

        let separated2 = separated_strings[1].clone();
        assert_eq!(
            separated2,
            line_string![
                (x: 0., y:1.0),
                (x: 0., y:2.0),
                (x: 0., y:3.0),
                (x: 0., y:4.0),
            ]
        );

        let separated3 = separated_strings[2].clone();
        assert_eq!(
            separated3,
            line_string![
                (x: 0., y:1.0),
                (x: 0., y:2.0),
                (x: 0., y:3.0),
            ]
        );

        let separated4 = separated_strings[3].clone();
        assert_eq!(
            separated4,
            line_string![
                (x: 0., y:1.0),
                (x: 0., y:2.0),
                (x: 0., y:3.0),
            ]
        );
    }
}
