//! Definition of the [`StationaryCurve`] type and the [`StationaryCurveSegment`] container type.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;

use colored::Colorize;
use contour::ContourBuilder;
use geo_types::{LineString, MultiLineString};
use nalgebra::Vector2;
use ndarray::{Array1, Array2, ArrayView1};
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use rsl_interpolation::Accelerator2d;

use dexter_machine::{
    EvalError, FluxCoordinateState, Machine, MachineError, MagneticFlux, MagneticFluxKind,
};

use crate::constants::{
    SC_CONTOUR_FLUX_POINTS, SC_CONTOUR_THETA_POINTS, SC_FLUX_BOUND, SC_THETA_BOUND,
    SC_THETA_CLAMP_RANGE, SCGD_FIRST_STEP, SCGD_MAX_DEVIATION, SCGD_MAX_ITERATIONS, SCGD_MAX_STEP,
    SCGD_MIN_STEP,
};

/// A container type that stores the θ and flux values of a single segment of the
/// [`StationaryCurve`].
///
/// The segment is firstly constructed by a contour generator, which is used as a set of initial
/// "guesses". Then, a gradient descent is performed at every point to reduce the deviation from
/// the true segment.
///
/// A gradient descent step is only accepted if the resulting deviation is less than the previous
/// one. Therefore, it is at least as accurate as the points calculated from the contour generator.
#[derive(Clone)]
pub struct StationaryCurveSegment {
    /// The `θ` array, in rads.
    theta: Array1<f64>,
    /// The flux array, in normalized units.
    flux: Array1<f64>,
    /// The kind of the magnetic flux w.r.t. which the segment is expressed.
    flux_kind: MagneticFluxKind,
    /// The value of the magnetic flux at the last closed surface.
    flux_last: f64,
    /// The average number of fitting iterations of every point.
    average_iterations: f64,
    /// The average number of successful descents (error reductions) of every point.
    average_descents: f64,
    /// The standard deviation of the segment at its construction from a contour line.
    start_std: f64,
    /// The standard deviation of the segment after the gradient descent.
    final_std: f64,
}

impl StationaryCurveSegment {
    /// Returns an [`ArrayView1`] over the `θ` values.
    #[must_use]
    pub fn flux(&self) -> ArrayView1<'_, f64> {
        self.flux.view()
    }

    /// Returns an [`ArrayView1`] over the flux values.
    #[must_use]
    pub fn theta(&self) -> ArrayView1<'_, f64> {
        self.theta.view()
    }

    /// Returns the length of the segment.
    #[must_use]
    #[expect(clippy::len_without_is_empty, reason = "segments cannot be empty")]
    pub fn len(&self) -> usize {
        self.theta.len()
    }

    /// Creates a `StationaryCurveSegment` from a [`geo_types::LineString`].
    ///
    /// The line string must be created from a contour generator with its origin at
    /// (`x0-dx/2`, `y0-dy/2`) and on a normalized grid. See [`StationaryCurve::build`].
    ///
    /// To correctly construct an open isoline, `line_string` must be a result of
    /// [`separate`].
    fn from_line_string(line_string: &LineString, flux_last: MagneticFlux) -> Self {
        let len = line_string.0.len();
        let mut theta = Vec::<f64>::with_capacity(len);
        let mut flux = Vec::<f64>::with_capacity(len);

        // NOTE: `coords()` are aligned correctly at the construction of the contour generator,
        // separated and clamped.
        line_string.coords().for_each(|coord| {
            theta.push(coord.x);
            flux.push(coord.y);
        });

        assert!(
            theta.iter().all(|theta| (0.0..TAU).contains(theta)),
            "Stationary Curve: Encountered 'θ' values outside the [0,2π] interval"
        );

        Self {
            flux_kind: flux_last.kind(),
            flux_last: flux_last.value(),
            theta: Array1::from_vec(theta),
            flux: Array1::from_vec(flux),
            average_iterations: 0.0,
            average_descents: 0.0,
            start_std: f64::NAN,
            final_std: f64::NAN,
        }
    }

    /// Calculates the standard deviation of the `𝜕B/𝜕𝜃` values evaluated on the segment's points.
    fn std(&self, machine: Machine) -> f64 {
        let acc = &mut Accelerator2d::new();
        let mut db = Vec::with_capacity(self.len());

        for index in 0..self.len() {
            let flux = self.flux_kind.to_magnetic_flux(self.flux[index]);
            let theta = self.theta[index];
            let Ok(value) = machine.bfield().eval_deriv_theta(flux, theta, acc) else {
                unreachable!("Points are in bounds by construction: (flux={flux:?}, theta={theta})")
            };
            db.push(value);
        }
        let db = Array1::from_vec(db);

        db.std(0.0)
    }

    /// Iterates over every point and performs a gradient descent.
    fn gradient_descent(&mut self, machine: Machine) {
        self.start_std = self.std(machine);
        for index in 0..self.len() {
            match self.gradient_descent_single_point(machine, index) {
                Ok(_) => (),
                Err(_) => unreachable!("Out of bounds evaluation in SC gradient descent"),
            };
        }
        self.final_std = self.std(machine);

        if self.start_std < self.final_std {
            let msg =
                "Stationary curve segment gradient descent made things worse".color("#ff9900");
            eprintln!("{msg}");
        }
    }

    /// Performs a gradient descent on every point of the segment.
    fn gradient_descent_single_point(
        &mut self,
        machine: Machine,
        index: usize,
    ) -> Result<(), EvalError> {
        let bfield = machine.bfield();

        let min_step = Vector2::new(SCGD_MIN_STEP.0 * TAU, SCGD_MIN_STEP.1 * self.flux_last);
        let max_step = Vector2::new(SCGD_MAX_STEP.0 * TAU, SCGD_MAX_STEP.1 * self.flux_last);
        let mut step = Vector2::new(SCGD_FIRST_STEP.0 * TAU, SCGD_FIRST_STEP.1 * self.flux_last);

        let acc = &mut Accelerator2d::new();
        let mut iterations = 0;
        let mut descents = 0;
        // TODO: This loop unnecessarily re-evaluates `bfield` quantities when the step was not accepted
        while iterations < SCGD_MAX_ITERATIONS {
            let old_theta = self.theta[index];
            let old_flux = self.flux_kind.to_magnetic_flux(self.flux[index]);

            let old_db = bfield.eval_deriv_theta(old_flux, old_theta, acc)?;
            let db_dtheta2 = bfield.eval_deriv_theta2(old_flux, old_theta, acc)?;
            let db_dflux2 = bfield.eval_deriv_flux2(old_flux, old_theta, acc)?;

            // Point towards the closest root of `𝜕B/𝜕𝜃`
            let grad_norm = Vector2::new(db_dtheta2, db_dflux2).normalize();
            let direction = if old_db.is_sign_positive() {
                -grad_norm
            } else {
                grad_norm
            };

            // Take a step towards the `𝜕B/𝜕𝜃 = 0` curve
            let new_theta = old_theta + step.x * direction.x;
            let new_flux = old_flux + step.y * direction.y;

            // Clamp the values within the bounds of the machine.
            // For `theta`, this means we avoid possible duplicate segments.
            // For `flux`, this means that we avoid getting too close to the axis and we prevent a
            // point from escaping the walls.
            let new_theta = new_theta.clamp(0.0, SC_THETA_BOUND.end);
            let new_flux = self.flux_kind.to_magnetic_flux(new_flux.value().clamp(
                SC_FLUX_BOUND.start * self.flux_last,
                SC_FLUX_BOUND.end * self.flux_last,
            ));

            let new_db = bfield.eval_deriv_theta(new_flux, new_theta, acc)?;

            // If the new point's `𝜕B/𝜕𝜃` evaluates closer to 0, accept it and continue with a
            // slightly smaller step. This guarantees that the result of the gradient descent is
            // at least as accurate as the initial contour line.
            //
            // If the new point lies further than the first one, discard it and continue with an
            // even smaller step.
            if new_db.abs() < old_db.abs() {
                self.theta[index] = new_theta;
                self.flux[index] = new_flux.value();
                step *= 0.9;
                step.x = step.x.min(max_step.x);
                step.y = step.y.min(max_step.y);
                descents += 1;
            } else {
                step *= 0.5;
                step.x = step.x.max(min_step.x);
                step.y = step.y.max(min_step.y);
            }

            // Break if we reached the desired accuracy.
            if new_db.abs() < SCGD_MAX_DEVIATION {
                break;
            }

            iterations += 1;
        }

        self.average_iterations += iterations as f64 / self.len() as f64;
        self.average_descents += descents as f64 / self.len() as f64;

        Ok(())
    }
}

/// The stationary curve of the Hamiltonian.
///
/// The stationary curve is defined through the equation `𝜕H/𝜕𝜃 = 0`. In the absence of an electric
/// field, this is equivalent to `𝜕B/𝜕𝜃 = 0`.
///
/// The stationary curve is described by one or more segments of the form `f(ψ,θ)=0`.
#[derive(Clone)]
#[non_exhaustive]
pub struct StationaryCurve {
    /// The kind of the magnetic flux w.r.t. which the curve is expressed.
    pub flux_kind: MagneticFluxKind,
    /// The curve's distinct segments.
    pub segments: Vec<StationaryCurveSegment>,
}

impl StationaryCurve {
    /// Builds a `StationaryCurve` from a [`Machine`].
    ///
    /// In the absence of an electric field, the stationary curve is built by calculating the
    /// separate `𝜕B/𝜕𝜃 = 0` segments on a contour of `𝜕B/𝜕𝜃`.
    ///
    /// # Errors
    ///
    /// No errors can occur from the function, but might do with the addition of the electric field.
    ///
    /// # Panics
    ///
    /// This function panics if [`ContourBuilder::lines`] returns an error.
    #[expect(clippy::panic_in_result_fn, reason = "unusual error, should be fatal")]
    pub fn build(machine: Machine) -> Result<Self, MachineError> {
        let (flux_kind, flux_last) = match machine.bfield().psi_state() {
            FluxCoordinateState::Good => (MagneticFluxKind::Toroidal, machine.qfactor().psi_last()),
            _ => (MagneticFluxKind::Poloidal, machine.qfactor().psip_last()),
        };
        let (theta_array, flux_array, grid) = Self::build_grid(machine);

        // ===========================

        // `dx` and `dy` defined as normal for an orthogonal grid
        let dx = theta_array[1] - theta_array[0];
        let dy = flux_array[1] - flux_array[0];

        // Origins must be shifted by half a step due to the way `ContourBuilder` defines the center
        // of the grid's squares
        let x0 = theta_array[0] - dx / 2.0;
        let y0 = flux_array[0] - dy / 2.0;

        let cb = ContourBuilder::new(theta_array.len(), flux_array.len(), true)
            .x_origin(x0)
            .y_origin(y0)
            .x_step(dx)
            .y_step(dy);

        let isolines = match cb.lines(grid.as_slice().expect("in logical order"), &[0.0]) {
            Ok(isolines) => isolines,
            Err(err) => {
                // If no contours are found, there is something wrong with the magnetic field
                panic!("ContourBuilder::lines panicked (error: {err})")
            }
        };

        // Flatten `lines: Vec<Line>` to a vector that holds references to the distinct
        // `LineStrings`, thus avoiding cloning. The line strings are then moved when
        // building the separateped strings.
        let strings: Vec<&LineString> = isolines
            .iter()
            .map(contour::Line::geometry)
            .collect::<Vec<&MultiLineString>>()
            .into_iter()
            .flatten()
            .collect();

        // Separate the discovered strings to individual lines.
        let mut separateped_strings: Vec<LineString> = strings
            .iter()
            .map(|string| separate(string, flux_array[flux_array.len() - 1]))
            .collect::<Vec<Vec<LineString>>>()
            .into_iter()
            .flatten()
            .collect();

        // Clamp the very slightly negative `thetas` that might occur due to the widened `theta` grid.
        separateped_strings.iter_mut().for_each(clamp_theta);

        // Create `StationaryCurveSegments` and perform the gradient descents
        let mut segments: Vec<StationaryCurveSegment> = separateped_strings
            .iter()
            .map(|string| StationaryCurveSegment::from_line_string(string, flux_last))
            .collect();

        segments
            .par_iter_mut()
            .for_each(|segment| segment.gradient_descent(machine));

        Ok(Self {
            flux_kind,
            segments,
        })
    }

    /// Builds the 2D grid on which the `StationaryCurve` is calculated.
    fn build_grid(machine: Machine) -> (Array1<f64>, Array1<f64>, Array2<f64>) {
        let (flux_kind, flux_last) = match machine.bfield().psi_state() {
            FluxCoordinateState::Good => (
                MagneticFluxKind::Toroidal,
                machine.qfactor().psi_last().value(),
            ),
            _ => (
                MagneticFluxKind::Poloidal,
                machine.qfactor().psip_last().value(),
            ),
        };

        // Widen the `theta` span to avoid edge phenomena
        let theta_first = -0.1;
        let theta_last = TAU + 0.1;
        // Avoid the singularity at the axis and getting too close to the last surface
        let flux_first = 1e-9 * flux_last;
        let flux_last = (1.0 - 1e-9) * flux_last;

        let theta_array = Array1::linspace(theta_first, theta_last, SC_CONTOUR_THETA_POINTS);
        let flux_array = Array1::linspace(flux_first, flux_last, SC_CONTOUR_FLUX_POINTS);
        let mut grid = Array2::<f64>::uninit((flux_array.len(), theta_array.len()));

        assert_eq!(grid.ncols(), theta_array.len(), "sanity check");
        assert_eq!(grid.nrows(), flux_array.len(), "sanity check");

        let acc = &mut Accelerator2d::new();

        for m in 0..grid.nrows() {
            let flux = flux_kind.to_magnetic_flux(flux_array[m]);
            for n in 0..grid.ncols() {
                let theta = theta_array[n].rem_euclid(TAU);
                if let Ok(db_dtheta) = machine.bfield().eval_deriv_theta(flux, theta, acc) {
                    grid[[m, n]] = MaybeUninit::new(db_dtheta)
                } else {
                    unreachable!("arrays are always in-bounds and the flux is always good")
                }
            }
        }

        // SAFETY: The loop passes from all elements and initializes them
        (theta_array, flux_array, unsafe { grid.assume_init() })
    }
}

/// Separates a [`LineString`] into multiple `LineStrings` by removing the parts that touch
/// the wall and lie outside the [0,2π] interval.
///
/// This is necessary as [`contour`] only yields closed surfaces. However, the "useless"
/// segments we want to discard always connect the true isolines by hugging the grid's
/// bounds (By construction of the contour lines through the marching squares algorigthm,
/// these "connecting" segments cannot lie anywhere but on the grid's bounds.). We can
/// use this to iterate over every string and safely separate these segments.
fn separate(string: &LineString, flux_last: f64) -> Vec<LineString> {
    let mut res = Vec::<LineString>::new();

    let xbound = SC_THETA_BOUND;
    let ybound = (SC_FLUX_BOUND.start * flux_last)..(SC_FLUX_BOUND.end * flux_last);

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
    res
}

/// Since we have widened the `theta` span, some points might land just below the `θ=0` line. We
/// clamp this values to `θ=0` exactly.
fn clamp_theta(string: &mut LineString) {
    string.coords_mut().for_each(|coord| {
        if SC_THETA_CLAMP_RANGE.contains(&coord.x) {
            coord.x = 0.0
        };
    });
}

impl std::fmt::Debug for StationaryCurveSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let last_index = self.len() - 1;
        let start = format_args!(
            "flux(normalized)={}, θ={}",
            self.flux[0] / self.flux_last,
            self.theta[0]
        );
        let end = format_args!(
            "flux(normalized)={}, θ={}",
            self.flux[last_index] / self.flux_last,
            self.theta[last_index]
        );
        let iterations = format_args!("{:.4}", self.average_iterations);
        let descents = format_args!("{:.4}", self.average_descents);
        f.debug_struct("StationaryCurveSegment")
            .field("start", &start)
            .field("end", &end)
            .field("length", &self.len())
            .field("average_iterations", &iterations)
            .field("average_descents", &descents)
            .field("start_std", &self.start_std)
            .field("final_std", &self.final_std)
            .finish()
    }
}

impl std::fmt::Debug for StationaryCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StationaryCurve")
            .field("flux_kind", &self.flux_kind)
            .field("number of segments", &self.segments.len())
            .field("segments", &self.segments)
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
