//! Defines the [`StationaryCurve`] type and its associated methods.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;
use std::ops::Range;
use std::time::{Duration, Instant};

use contour::{ContourBuilder, Line};
use dexter_machine::{Accelerator2d, Machine, MagneticFluxArray, MagneticFluxKind};
use geo_types::{LineString, MultiLineString};
use ndarray::{Array1, Array2, ArrayRef1};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{CriticalPoint, ParametrizedSegment};

/// The density of the (θ, ψ/ψp) grid on which to search the stationary curve contours.
const SC_CONTOUR_DENSITY: (usize, usize) = (300, 300);

// ===============================================================================================

/// The stationary curve of the Hamiltonian.
///
/// The stationary curve is defined through the equation `𝜕H/𝜕𝜃 = 0`. In the absence of an electric
/// field, this is equivalent to `𝜕B/𝜕𝜃 = 0`.
///
/// The stationary curve is described by one or more [`parametrized segments`][ParametrizedSegment]
/// of the form `𝜕H/𝜕𝜃(ψ,θ)=0`.
///
/// This type is built through
/// [`GCHamiltonian::build_stationary_curve`][crate::GCHamiltonian::build_stationary_curve].
#[derive(Clone)]
pub struct StationaryCurve {
    /// The separated segments comprising the `StationaryCurve`.
    segments: Box<[ParametrizedSegment]>,
    /// The duration of the [`Self::calculate_critical_points`] routine.
    duration: Option<Duration>,
}

impl StationaryCurve {
    /// Returns a reference to the stored [`ParametrizedSegments`][ParametrizedSegment].
    #[must_use]
    pub fn segments(&self) -> &[ParametrizedSegment] {
        &self.segments
    }

    /// Returns the number of stored [`ParametrizedSegments`][ParametrizedSegment].
    #[must_use]
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Returns `true` if the curve contains no [`ParametrizedSegments`][ParametrizedSegment].
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
    ///
    /// # Process
    ///
    /// 1. Build the (θ, ψ/ψp) *orthogonal* grid with [`build_grid`].
    /// 2. Set up the [`ContourBuilder`].
    /// 3. Extract [`ContourBuilder::lines`] into a [`Box<ParametrizedSegment>`] and store it to `self`.
    ///
    /// # Panics
    ///
    /// This method panics if [`ContourBuilder::lines`] returns an error. This is unusual since at
    /// least one segment of `𝜕B/𝜕𝜃 = 0` should always exist and therefore should be fatal.
    pub(crate) fn build(machine: Machine) -> Self {
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
        let segments = separated_strings_into_segments(separated_strings, machine.good_flux());

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
fn isolines_to_strings(isolines: &[Line]) -> Box<[&LineString]> {
    isolines
        .iter()
        .map(Line::geometry)
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
) -> Box<[ParametrizedSegment]> {
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

            ParametrizedSegment::new(flux_array, theta_array)
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

impl std::fmt::Debug for StationaryCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StationaryCurve")
            .field(
                "flux_kind",
                &self.segments.first().map(ParametrizedSegment::flux_kind),
            )
            .field("number of segments", &self.segments.len())
            .field("segments", &self.segments)
            .field("duration", &self.duration)
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
