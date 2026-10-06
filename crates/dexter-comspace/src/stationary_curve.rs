//! Definition of the [`StationaryCurve`] type and the [`StationaryCurveSegment`] container type.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;
use std::ops::Range;

use contour::ContourBuilder;
use geo_types::{LineString, MultiLineString};
use ndarray::{Array1, Array2, ArrayView1};
use rsl_interpolation::Accelerator2d;

use crate::constants::{SC_CONTOUR_FLUX_POINTS, SC_CONTOUR_THETA_POINTS};
use dexter_machine::{FluxCoordinateState, Machine, MachineError, MagneticFluxKind};

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

    /// Creates a `StationaryCurveSegment` from a [`geo_types::LineString`].
    ///
    /// The line string must be created from a contour generator with its origin at (0.0) and its
    /// dimensions normalized to 2π and the last closed surface accordingly.
    ///
    /// To correctly construct an open isoline, `line_string` must be a result of
    /// [`StationaryCurve::separate`].
    fn from_line_string(line_string: &LineString) -> Self {
        let len = line_string.0.len();
        let mut theta = Vec::<f64>::with_capacity(len);
        let mut flux = Vec::<f64>::with_capacity(len);

        line_string.coords().for_each(|coord| {
            theta.push(coord.x);
            flux.push(coord.y);
        });

        Self {
            theta: Array1::from_vec(theta),
            flux: Array1::from_vec(flux),
        }
    }
}

/// The stationary curve of a Hamiltonian.
///
/// The stationary curve is defined through the equation `𝜕H/𝜕𝜃 = 0`. In the absence of an electric
/// field, this is equivalent to `𝜕B/𝜕𝜃 = 0`.
///
/// The stationary curve is described by one or more segments of the form `f(ψ,θ)=0`.
#[derive(Debug, Clone)]
pub struct StationaryCurve {
    /// The curve's `μ` parameter.
    ///
    /// This value is currently unused, but will be utilized when the electric field is added.
    pub _mu: f64,
    /// The curve's `Pζ` parameter.
    ///
    /// This value is currently unused, but will be utilized when the electric field is added.
    pub _pzeta: f64,
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
        let (flux_kind, flux_last) = match machine.qfactor().psi_state() {
            FluxCoordinateState::Good => (
                MagneticFluxKind::Toroidal,
                machine.qfactor().psi_last().value(),
            ),
            _ => (
                MagneticFluxKind::Poloidal,
                machine.qfactor().psip_last().value(),
            ),
        };
        let (theta_array, flux_array, grid) = Self::build_grid(machine);

        let flux_bound = flux_array.last().expect("non-empty by construction");
        let bounds = (1e-10 * flux_bound)..((1.0 - 1e-10) * flux_bound);

        // ===========================

        let cb = ContourBuilder::new(theta_array.len(), flux_array.len(), true)
            .x_origin(0.0)
            .y_origin(0.0)
            .x_step(TAU / SC_CONTOUR_THETA_POINTS as f64)
            .y_step(flux_last / SC_CONTOUR_FLUX_POINTS as f64);

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

        let separateped_strings: Vec<LineString> = strings
            .iter()
            .map(|string| Self::separate(string, &bounds))
            .collect::<Vec<Vec<LineString>>>()
            .into_iter()
            .flatten()
            .collect();

        let segments: Vec<StationaryCurveSegment> = separateped_strings
            .iter()
            .map(StationaryCurveSegment::from_line_string)
            .collect();

        Ok(Self {
            flux_kind,
            segments,
            _mu: f64::NAN,
            _pzeta: f64::NAN,
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

        // NOTE: `Array1::linspace` includes values from `start` to `end` *inclusive*.
        // We prevent the `θ` array from reaching 2π as that might yield duplicate segments.
        // We prevent the flux array from going too close to the axis as it may yield multiple
        // length-1 duplicate segments. However, we allow it to reach the last closed flux surface
        // so we can safely separate the segments later.
        let flux_array = Array1::linspace(1e-10 * flux_last, flux_last, SC_CONTOUR_FLUX_POINTS);
        let theta_array = Array1::linspace(0.0, (1.0 - 1e-10) * TAU, SC_CONTOUR_THETA_POINTS);
        let mut grid = Array2::<f64>::uninit((flux_array.len(), theta_array.len()));

        assert_eq!(grid.ncols(), theta_array.len(), "sanity check");
        assert_eq!(grid.nrows(), flux_array.len(), "sanity check");

        let acc = &mut Accelerator2d::new();

        for m in 0..grid.nrows() {
            let flux = flux_kind.to_magnetic_flux(flux_array[m]);
            for n in 0..grid.ncols() {
                let theta = theta_array[n];
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

    /// Separates a [`LineString`] into multiple `LineStrings` by removing the parts that touch
    /// the wall.
    ///
    /// This is necessary as [`contour`] only yields closed surfaces. However, the "useless"
    /// segments we want to discard always connect the true isolines by hugging the grid's
    /// bounds (By construction of the contour lines through the marching squares algorigthm,
    /// these "connecting" segments cannot lie anywhere but on the grid's bounds.). We can
    /// use this to iterate over every string and safely separate these segments.
    fn separate(string: &LineString, bounds: &Range<f64>) -> Vec<LineString> {
        let mut res = Vec::<LineString>::new();

        // Keeps track of the end index of every discovered valid segment
        let mut end = 0;

        while end < string.0.len() {
            let iterator = string
                .coords()
                .skip(end)
                .take_while(|coord| {
                    end += 1;
                    bounds.contains(&coord.y)
                })
                .copied();
            let new_line = LineString::from_iter(iterator);
            if new_line.0.len() > 1 {
                res.push(new_line);
            }
        }

        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use geo_types::line_string;

    const BOUNDS: Range<f64> = 0.1..10.0;

    fn assert_separated_equals_expected(string: LineString, expected: LineString) {
        let separated_strings = StationaryCurve::separate(&string, &BOUNDS);
        assert_eq!(separated_strings.len(), 1);
        let separated = separated_strings[0].clone();
        assert_eq!(separated, expected);
    }

    #[test]
    fn separate_line_string_none1() {
        let string: LineString = line_string![
            (x: 0., y:0.),
            (x: 0., y:0.),
        ];
        let separateped_strings = StationaryCurve::separate(&string, &BOUNDS);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none2() {
        let string: LineString = line_string![
            (x: 0., y:20.0),
            (x: 0., y:30.0),
        ];
        let separateped_strings = StationaryCurve::separate(&string, &BOUNDS);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none3() {
        let string: LineString = line_string![
            (x: 0., y:3.0),
        ];
        let separateped_strings = StationaryCurve::separate(&string, &BOUNDS);
        assert_eq!(separateped_strings.len(), 0);
    }

    #[test]
    fn separate_line_string_none4() {
        let string: LineString = line_string![
            (x: 0., y:3.0),
            (x: 0., y:0.),
        ];
        let separateped_strings = StationaryCurve::separate(&string, &BOUNDS);
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
        let separated_strings = StationaryCurve::separate(&string, &BOUNDS);
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
            (x: 0., y:100.0),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
            (x: 0., y:4.0),
            (x: 0., y:100.0),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
            (x: 0., y:0.),
            (x: 0., y:1.0),
            (x: 0., y:2.0),
            (x: 0., y:3.0),
        ];
        let separated_strings = StationaryCurve::separate(&string, &BOUNDS);
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
