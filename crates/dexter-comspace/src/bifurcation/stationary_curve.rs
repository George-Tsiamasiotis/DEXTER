//! Definition of the [`StationaryCurve`] type and the [`StationaryCurveSegment`] container type.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;
use std::ops::Range;

use approx::relative_eq;
use contour::ContourBuilder;
use geo_types::{LineString, MultiLineString};
use ndarray::{Array1, Array2, ArrayView1};
use rsl_interpolation::Accelerator2d;

use crate::constants::{SC_CONTOUR_FLUX_POINTS, SC_CONTOUR_THETA_POINTS};
use dexter_machine::{FluxCoordinateState, Machine, MachineError, MagneticFluxKind};

/// The bounds of the `theta` array in which a stationary curve segment is considered valid.
const THETA_BOUND: Range<f64> = (-1e-4 * TAU)..((1.0 - 1e-4) * TAU);

/// The bounds of the `flux` array in which a stationary curve segment is considered valid.
/// Units are normalized w.r.t. the flux value at the last closed surface.
const FLUX_BOUND: Range<f64> = (1e-8)..(1.0 - 1e-10);

/// The `theta` range in which slightly negative values are clamped to 0.
const THETA_CLAMP_RANGE: Range<f64> = -1e-4 * TAU..0.0;

/// A container type that stores the θ and flux values of a single segment of the stationary curve.
#[derive(Debug, Clone)]
pub struct StationaryCurveSegment {
    /// The `θ` array, in rads.
    theta: Array1<f64>,
    /// The flux array, in normalized units.
    flux: Array1<f64>,
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
    /// [`StationaryCurve::separate`].
    fn from_line_string(line_string: &LineString) -> Self {
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
            theta: Array1::from_vec(theta),
            flux: Array1::from_vec(flux),
        }
    }

    /// Checks if the calculated `theta` and `flux` arrays indeed evaluate to `𝜕B/𝜕𝜃 = 0`,
    /// panicking otherwise.
    ///
    /// Note that the tolerance is set relatively low (1%) due to inaccuracies close to the wall and
    /// the magnetic axis.
    fn check_validity(&self, machine: Machine, flux_kind: &MagneticFluxKind) {
        let acc = &mut Accelerator2d::new();
        for index in 0..self.len() {
            let flux = flux_kind.to_magnetic_flux(self.flux[index]);
            let theta = self.theta[index];
            let Ok(db) = machine.bfield().eval_deriv_theta(flux, theta, acc) else {
                unreachable!("segments are always in bounds")
            };

            // TODO: make stationary curve more accurate
            assert!(
                relative_eq!(db, 0.0, epsilon = 1e-2),
                "Stationary curve: 𝜕B/𝜕𝜃(flux={flux:?}, theta={theta}) = {db}"
            )
        }
    }
}

/// The stationary curve of the Hamiltonian.
///
/// The stationary curve is defined through the equation `𝜕H/𝜕𝜃 = 0`. In the absence of an electric
/// field, this is equivalent to `𝜕B/𝜕𝜃 = 0`.
///
/// The stationary curve is described by one or more segments of the form `f(ψ,θ)=0`.
#[derive(Debug, Clone)]
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
        let flux_kind = match machine.bfield().psi_state() {
            FluxCoordinateState::Good => MagneticFluxKind::Toroidal,
            _ => MagneticFluxKind::Poloidal,
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

        // Create `StationaryCurveSegments` and check them
        let segments: Vec<StationaryCurveSegment> = separateped_strings
            .iter()
            .map(StationaryCurveSegment::from_line_string)
            .collect();
        segments
            .iter()
            .for_each(|segment| segment.check_validity(machine, &flux_kind));

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

    let xbound = THETA_BOUND;
    let ybound = (FLUX_BOUND.start * flux_last)..(FLUX_BOUND.end * flux_last);

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
        if THETA_CLAMP_RANGE.contains(&coord.x) {
            coord.x = 0.0
        };
    });
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
