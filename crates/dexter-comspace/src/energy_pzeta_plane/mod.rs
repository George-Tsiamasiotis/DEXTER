//! Representation of the COM space `(E, Pζ, μ=const)`.

#![expect(clippy::min_ident_chars, reason = "parabola a, b, c coefficients")]

mod orbit_type;
mod tp_boundary;

use std::time::{Duration, Instant};

use dexter_machine::MagneticFlux::*;
use ndarray::{ArrayRef1, ArrayView1};
use parabola::Parabola;
use rsl_interpolation::Accelerator2d;

use dexter_machine::Machine;

pub use orbit_type::{EnergyPzetaPosition, OrbitType};
pub use tp_boundary::TrappedPassingBoundary;

use crate::{Boundary, BoundaryExtremum, CriticalPoint, StationaryCurve};

/// Representation of the COM space `(E, Pζ, μ=const)`.
///
/// ![Parabolas][parabolas].
#[cfg_attr(
    feature="doc-images",
    doc = ::embed_doc_image::embed_image!("parabolas", "../../docs/assets/parabolas.svg"))
]
#[derive(Clone)]
pub struct EnergyPzetaPlane {
    /// The constant magnetic moment `μ`.
    mu: f64,
    /// The magnetic axis (MA) parabola.
    axis_parabola: Parabola,
    /// The left wall (LW) parabola.
    left_wall_parabola: Parabola,
    /// The right wall (RW) parabola.
    right_wall_parabola: Parabola,
    /// The trapped-passing boundary curves.
    tp_boundary: TrappedPassingBoundary,
    /// The calculated stationary curve.
    ///
    /// Even though the stationary curve does not depend on mu, it is handy it to store it here. It
    /// is only calculated when [`Self::calculate_critical_points`] is called.
    stationary_curve: Option<StationaryCurve>,
    /// The calculated boundary of the Hamiltonian.
    boundary: Option<Boundary>,
    /// The calculated critical points.
    critical_points: Box<[CriticalPoint]>,
    /// The calculated boundary extremum points.
    boundary_extrema: Box<[BoundaryExtremum]>,
    /// The duration of the [`Self::calculate_critical_points`] subroutine.
    critical_point_search_duration: Option<Duration>,
    /// The duration of the [`Self::calculate_boundary_extrema`] subroutine.
    boundary_extremum_search_duration: Option<Duration>,
    /// The duration of the [`Self::classify_critical_points`] subroutine.
    critical_point_classification_duration: Option<Duration>,
    /// The duration of the [`Self::classify_boundary_extrema`] subroutine.
    boundary_extremum_classification_duration: Option<Duration>,
}

impl EnergyPzetaPlane {
    /// Creates a new `EnergyPzetaPlane` from a set magnetic moment `μ=const` value.
    #[must_use]
    pub fn from_mu(machine: Machine, mu: f64) -> Self {
        Self {
            mu,
            axis_parabola: Self::build_magnetic_axis_parabola(machine, mu),
            left_wall_parabola: Self::build_left_wall_parabola(machine, mu),
            right_wall_parabola: Self::build_right_wall_parabola(machine, mu),
            tp_boundary: TrappedPassingBoundary::new(machine, mu),
            stationary_curve: None,
            boundary: None,
            critical_points: Box::new([]),
            boundary_extrema: Box::new([]),
            critical_point_search_duration: None,
            boundary_extremum_search_duration: None,
            critical_point_classification_duration: None,
            boundary_extremum_classification_duration: None,
        }
    }
}

/// Bifurcation Analysis.
impl EnergyPzetaPlane {
    /// Performs a full Bifurcation Analysis on the `E-Pζ` plane.
    ///
    /// This is a convenience method that simply calls the following methods:
    ///
    /// + [`Self::calculate_critical_points`],
    /// + [`Self::calculate_boundary_extrema`],
    /// + [`Self::classify_critical_points`],
    /// + [`Self::classify_boundary_extrema`].
    pub fn bifurcation_analysis(&mut self, machine: Machine, pzeta_array: &ArrayRef1<f64>) {
        self.calculate_critical_points(machine, pzeta_array);
        self.calculate_boundary_extrema(machine, pzeta_array);
        self.classify_critical_points(machine);
        self.classify_boundary_extrema(machine);
    }

    /// Calculates the [`critical points`][CriticalPoint] of the `(E, Pζ; μ)` slice over all values
    /// of `pzeta_array` and stores them.
    pub fn calculate_critical_points(&mut self, machine: Machine, pzeta_array: &ArrayRef1<f64>) {
        let start = Instant::now();
        let mut curve = StationaryCurve::new(machine);
        self.critical_points = curve.calculate_critical_points(machine, pzeta_array, self.mu);
        self.critical_point_search_duration = Some(start.elapsed());
        self.stationary_curve = Some(curve);
    }

    /// Calculates the [`boundary_extrema`][BoundaryExtremum] of the `(E, Pζ; μ)` slice over all
    /// values of `pzeta_array` and stores them.
    pub fn calculate_boundary_extrema(&mut self, machine: Machine, pzeta_array: &ArrayRef1<f64>) {
        let start = Instant::now();
        let boundary = Boundary::new(machine);
        self.boundary_extrema = boundary.calculate_boundary_extrema(machine, pzeta_array, self.mu);
        self.boundary_extremum_search_duration = Some(start.elapsed());
        self.boundary = Some(boundary);
    }

    /// Classifies the already calculated critical points into one of
    /// [`CriticalPointKind`][crate::CriticalPointKind].
    pub fn classify_critical_points(&mut self, machine: Machine) {
        let start = Instant::now();
        self.critical_points
            .iter_mut()
            .for_each(|point| point.classify(machine));
        self.critical_point_classification_duration = Some(start.elapsed());
    }

    /// Classifies the already calculated boundary extrema into one of
    /// [`BoundaryExtremumKind`][crate::BoundaryExtremumKind].
    pub fn classify_boundary_extrema(&mut self, machine: Machine) {
        let start = Instant::now();
        self.boundary_extrema
            .iter_mut()
            .for_each(|point| point.classify(machine));
        self.boundary_extremum_classification_duration = Some(start.elapsed());
    }
}

// Parabola builders
impl EnergyPzetaPlane {
    /// Constructs the Magnetic Axis (MA) parabola:
    ///
    /// E(Pζ) = Pζ²B²/2g² + μB + Φ.
    ///
    /// where all values are evaluated at `(ψ, θ) = (0, 0)`.
    #[must_use]
    fn build_magnetic_axis_parabola(machine: Machine, mu: f64) -> Parabola {
        let acc = &mut Accelerator2d::new();
        let psit_axis = Toroidal(0.0);
        let psip_axis = Poloidal(0.0);

        // Use `unwrap_or_else` for lazy evaluation and avoid panic
        let gaxis = machine
            .current()
            .eval_g(psit_axis, acc.xacc())
            .unwrap_or_else(|_| {
                machine
                    .current()
                    .eval_g(psip_axis, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let baxis = machine
            .bfield() // This might be redundant
            .eval_b(psit_axis, 0.0, acc)
            .unwrap_or_else(|_| {
                machine
                    .bfield()
                    .eval_b(psip_axis, 0.0, acc)
                    .expect("At least one of the evaluations will always succeed")
            });
        let phi_axis = machine.efield().eval_phi(psit_axis).unwrap_or_else(|_| {
            machine
                .efield()
                .eval_phi(psip_axis)
                .expect("At least one of the evaluations will always succeed")
        });

        Parabola {
            a: (baxis / gaxis).powi(2) / 2.0,
            b: 0.0,
            c: mu * baxis + phi_axis,
        }
    }

    /// Constructs the Left Wall (LW) parabola:
    ///
    /// E(Pζ) = (Pζ + ψp)²B²/2g² + μB + Φ.
    ///
    /// where all values are evaluated at `(ψ, θ) = (ψlast, π)`.
    #[must_use]
    fn build_left_wall_parabola(machine: Machine, mu: f64) -> Parabola {
        use std::f64::consts::PI;
        let psit_last = machine.qfactor().psit_last();
        let psip_last = machine.qfactor().psip_last();
        let acc = &mut Accelerator2d::new();

        // Use `unwrap_or_else` for lazy evaluation and avoid panic
        let glast = machine
            .current()
            .eval_g(psit_last, acc.xacc())
            .unwrap_or_else(|_| {
                machine
                    .current()
                    .eval_g(psip_last, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = machine
            .bfield()
            .eval_b(psit_last, PI, acc)
            .unwrap_or_else(|_| {
                machine
                    .bfield()
                    .eval_b(psip_last, PI, acc)
                    .expect("At least one of the evaluations will always succeed")
            });
        let phi_last = machine.efield().eval_phi(psit_last).unwrap_or_else(|_| {
            machine
                .efield()
                .eval_phi(psip_last)
                .expect("At least one of the evaluations will always succeed")
        });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last.value();
        let k = mu * blast + phi_last;
        Parabola::from_square(a, h, k)
    }

    /// Constructs the Right Wall (RW) parabola:
    ///
    /// E(Pζ) = (Pζ + ψp)²B²/2g² + μB + Φ.
    ///
    /// where all values are evaluated at `(ψ, θ) = (ψlast, 0)`.
    #[must_use]
    fn build_right_wall_parabola(machine: Machine, mu: f64) -> Parabola {
        let psit_last = machine.qfactor().psit_last();
        let psip_last = machine.qfactor().psip_last();
        let acc = &mut Accelerator2d::new();

        // Use `unwrap_or_else` for lazy evaluation and avoid panic
        let glast = machine
            .current()
            .eval_g(psit_last, acc.xacc())
            .unwrap_or_else(|_| {
                machine
                    .current()
                    .eval_g(psip_last, acc.xacc())
                    .expect("At least one of the evaluations will always succeed")
            });
        let blast = machine
            .bfield()
            .eval_b(psit_last, 0.0, acc)
            .unwrap_or_else(|_| {
                machine
                    .bfield()
                    .eval_b(psip_last, 0.0, acc)
                    .expect("At least one of the evaluations will always succeed")
            });
        let phi_last = machine.efield().eval_phi(psit_last).unwrap_or_else(|_| {
            machine
                .efield()
                .eval_phi(psip_last)
                .expect("At least one of the evaluations will always succeed")
        });

        let a = (blast / glast).powi(2) / 2.0;
        let h = psip_last.value();
        let k = mu * blast + phi_last;
        Parabola::from_square(a, h, k)
    }
}

impl EnergyPzetaPlane {
    /// Returns a reference to the magnetic axis (MA) [`Parabola`].
    #[must_use]
    pub fn axis_parabola(&self) -> &Parabola {
        &self.axis_parabola
    }

    /// Returns a reference to the left wall (LW) [`Parabola`].
    #[must_use]
    pub fn left_wall_parabola(&self) -> &Parabola {
        &self.left_wall_parabola
    }

    /// Returns a reference to the right wall (RW) [`Parabola`].
    #[must_use]
    pub fn right_wall_parabola(&self) -> &Parabola {
        &self.right_wall_parabola
    }

    /// Returns a reference to the [`TrappedPassingBoundary`].
    #[must_use]
    pub fn tp_boundary(&self) -> &TrappedPassingBoundary {
        &self.tp_boundary
    }

    /// Returns the constant magnetic moment `μ`.
    #[must_use]
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// Returns the [`TrappedPassingBoundary`]'s `Pζ = [-ψp_last, 0]` interval array.
    #[must_use]
    pub fn tp_pzeta_values(&self) -> ArrayView1<'_, f64> {
        self.tp_boundary.pzeta()
    }

    /// Returns the [`TrappedPassingBoundary`]'s upper curve.
    #[must_use]
    pub fn tp_upper_values(&self) -> ArrayView1<'_, f64> {
        self.tp_boundary.upper()
    }

    /// Returns the [`TrappedPassingBoundary`]'s lower curve.
    #[must_use]
    pub fn tp_lower_values(&self) -> ArrayView1<'_, f64> {
        self.tp_boundary.lower()
    }

    /// Returns a reference on the stored [`CriticalPoints`][CriticalPoint].
    #[must_use]
    pub fn critical_points(&self) -> &[CriticalPoint] {
        &self.critical_points
    }

    /// Returns a reference on the stored [`BoundaryExtrema`][BoundaryExtremum].
    #[must_use]
    pub fn boundary_extrema(&self) -> &[BoundaryExtremum] {
        &self.boundary_extrema
    }

    /// Returns an array with the calculated critical points' `Pζ`.
    ///
    /// The points are visited in order of creation, as in [`Self::critical_points`].
    #[must_use]
    pub fn critical_points_pzetas(&self) -> Box<[f64]> {
        self.critical_points
            .iter()
            .map(CriticalPoint::pzeta)
            .collect()
    }

    /// Returns an array with the calculated critical points' energies.
    ///
    /// The points are visited in order of creation, as in [`Self::critical_points`].
    #[must_use]
    pub fn critical_points_energies(&self) -> Box<[f64]> {
        self.critical_points
            .iter()
            .map(CriticalPoint::energy)
            .collect()
    }

    /// Returns an array with the calculated boundary extrema's `Pζ`.
    ///
    /// The points are visited in order of creation, as in [`Self::boundary_extrema`].
    #[must_use]
    pub fn boundary_extrema_pzetas(&self) -> Box<[f64]> {
        self.boundary_extrema
            .iter()
            .map(BoundaryExtremum::pzeta)
            .collect()
    }

    /// Returns an array with the calculated boundary extrema's energies.
    ///
    /// The points are visited in order of creation, as in [`Self::boundary_extrema`].
    #[must_use]
    pub fn boundary_extrema_energies(&self) -> Box<[f64]> {
        self.boundary_extrema
            .iter()
            .map(BoundaryExtremum::energy)
            .collect()
    }

    /// Prints the durations of the Bifurcation analysis' subroutines.
    pub fn print_durations(&self) {
        println!("Critical points: {}", self.critical_points.len());
        println!("Boundary extrema: {}", self.boundary_extrema.len());
        if let Some(str) = self.critical_point_search_duration {
            println!("Critical point search duration: {str:?}")
        };
        if let Some(str) = self.boundary_extremum_search_duration {
            println!("Boundary extremum search duration: {str:?}")
        };
        if let Some(str) = self.critical_point_classification_duration {
            println!("Critical point classification duration: {str:?}")
        };
        if let Some(str) = self.boundary_extremum_classification_duration {
            println!("Boundary extremum classification duration: {str:?}")
        };
    }
}

impl std::fmt::Debug for EnergyPzetaPlane {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EnergyPzetaPlane")
            .field("mu", &self.mu)
            .field("axis_parabola", &self.axis_parabola)
            .field("left_wall_parabola", &self.left_wall_parabola)
            .field("right_wall_parabola", &self.right_wall_parabola)
            .field("stationary_curve", &self.stationary_curve)
            .field("critical_points", &self.critical_points.len())
            .field("boundary_extrema", &self.boundary_extrema.len())
            .finish()
    }
}
