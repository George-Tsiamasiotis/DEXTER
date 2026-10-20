//! Defines the [`Boundary`] type and the methods for calculating its [`BoundaryExtremum`] points.

use std::f64::consts::TAU;
use std::mem::MaybeUninit;
use std::time::{Duration, Instant};

use dexter_machine::{Accelerator, Accelerator2d, Machine};
use ndarray::{Array1, ArrayRef1};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use roots::{Convergency, find_root_brent};
use rsl_interpolation::{Interpolation, LinearInterpolator};

use crate::{BoundaryExtremum, gc_hamiltonian::ParametrizedSegment};

// ===============================================================================================

/// Convergency criteria for finding the root of `𝜕B/𝜕𝜃`.
struct BoundaryExtremumConvergency;

impl BoundaryExtremumConvergency {
    /// Maximum allowed value for `𝜕B/𝜕𝜃`.
    const EPS: f64 = 1e-16;
    /// Maximum allowed iterations.
    const MAX_ITER: usize = 50;
    /// Maximum allowed difference between the 2 last `t` values.
    const MAX_PARAMETER_DIFF: f64 = 1e-7;
}

#[expect(clippy::renamed_function_params, reason = "clarity")]
impl Convergency<f64> for BoundaryExtremumConvergency {
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

/// Representation of the Hamiltonian's boundary surface (e.g the Last Closed Flux Surface).
#[derive(Debug, Clone)]
#[doc(alias = "last_closed_flux_surface")]
pub struct Boundary {
    /// The (t, ψ, θ) parametrized segment.
    segment: ParametrizedSegment,
    /// The duration of the [`Self::calculate_boundary_extrema`] routine.
    duration: Option<Duration>,
}

impl Boundary {
    /// The length of the contained parametrized segment.
    const LENGTH: usize = 400;

    /// Creates a new `BoundarySegment`.
    #[must_use]
    pub(crate) fn new(machine: Machine) -> Self {
        let flux_last = machine.last(machine.good_flux());
        let flux_values = Array1::from_elem(Self::LENGTH, flux_last.value() * (1.0 - 1e-9));
        let flux_array = flux_last.kind().to_magnetic_flux_array(flux_values);
        let theta_array = Array1::linspace(0.0, TAU + 0.1, Self::LENGTH); // Include roots at 2π
        let segment = ParametrizedSegment::new(flux_array, theta_array);
        Self {
            segment,
            duration: None,
        }
    }

    /// Calculates the boundary extrema of the `(E, Pζ; μ)` slice over all values of `pzeta_array` and
    /// stores them.
    pub(crate) fn calculate_boundary_extrema(
        &mut self,
        machine: Machine,
        pzeta_array: &ArrayRef1<f64>,
        mu: f64,
    ) -> Box<[BoundaryExtremum]> {
        let start = Instant::now();

        let db_dtheta_array = self.calculate_db_dtheta(machine);
        let root_indices = ParametrizedSegment::find_root_indices(&db_dtheta_array);
        let exact_roots = self.find_exact_roots(machine, &root_indices);

        let slice = pzeta_array.as_slice().unwrap_or(&[]); // do nothing
        let points: Vec<BoundaryExtremum> = slice
            .into_par_iter()
            .map(|pzeta| self.extrema_from_exact_roots(machine, &exact_roots, *pzeta, mu))
            .flatten()
            .collect();

        self.duration = Some(start.elapsed());
        points.into_boxed_slice()
    }

    /// Calculates the `𝜕B/𝜕𝜃` from the contained `flux` and `theta` arrays and the constants
    /// of motion `Pζ` and `μ`.
    #[must_use]
    fn calculate_db_dtheta(&self, machine: Machine) -> Array1<f64> {
        let acc = &mut Accelerator2d::new();
        let Some(flux) = self.segment.flux_array.get(0) else {
            unreachable!("segment cannot be empty")
        };

        let mut array = Array1::<f64>::uninit(self.segment.len());

        for n in 0..array.len() {
            let theta = self.segment.theta_array[n].rem_euclid(TAU);
            if let Ok(db_dtheta) = machine.bfield().eval_deriv_theta(flux, theta, acc) {
                array[n] = MaybeUninit::new(db_dtheta)
            } else {
                unreachable!("arrays are in-bounds and flux is Good")
            }
        }

        // SAFETY: The loop passes from all elements and initializes them
        unsafe { array.assume_init() }
    }

    /// Finds the true roots of `𝜕B/𝜕𝜃`.
    ///
    /// For every `t[i]`, `t[i+1]` interval where `𝜕B/𝜕𝜃` changes sign, we construct the 1D function:
    /// `f(t) = f(ψLCFS, θ(t)) = 𝜕B/𝜕𝜃(t)` and find its root.
    pub(crate) fn find_exact_roots(&self, machine: Machine, root_indices: &[usize]) -> Box<[f64]> {
        let mut exact_roots = Vec::new();
        root_indices.iter().for_each(|index| {
            // Index into the arrays at the edges of the `[i, i+1]` interval. By construction of
            // the parameter `t`, it always holds that `tlo<thi` and therefore the interpolation
            // is always defined.
            let tlo = self.segment.t[*index]; // == index
            let thi = self.segment.t[*index + 1]; // == index + 1
            let theta_lo = self.segment.theta_array[*index];
            let theta_hi = self.segment.theta_array[*index + 1];
            let ta = [tlo, thi];

            let ya_theta = [theta_lo, theta_hi];
            let Some(flux) = self.segment.flux_array.get(0) else {
                unreachable!("segment cannot be empty")
            };

            let linear_acc = &mut Accelerator::new();
            let b_acc = &mut Accelerator2d::new();

            let fun = |t: f64| {
                if !(tlo..=thi).contains(&t) {
                    return f64::MAX; // stop the solver. Landing here should be extremely unlikely
                }

                // `t -> θ` interpolation
                let theta = LinearInterpolator
                    .eval(&ta, &ya_theta, t, linear_acc)
                    .expect("arrays and values ok");

                machine
                    .bfield()
                    .eval_deriv_theta(flux, theta, b_acc)
                    .unwrap_or(f64::NAN) // stops the routine
            };

            let Ok(t_root) = find_root_brent(tlo, thi, fun, &mut BoundaryExtremumConvergency)
            else {
                return;
            };

            exact_roots.push(t_root);
        });
        exact_roots.into_boxed_slice()
    }

    /// Calculates the boundary extrema by interpolating each exact root to the corresponding `θ` and
    /// `ψ/ψp` points.
    pub(crate) fn extrema_from_exact_roots(
        &self,
        machine: Machine,
        exact_roots: &[f64],
        pzeta: f64,
        mu: f64,
    ) -> Vec<BoundaryExtremum> {
        self.segment
            .interpolate_exact_roots(exact_roots)
            .into_iter()
            .filter_map(|(theta, flux)| BoundaryExtremum::new(machine, theta, flux, pzeta, mu).ok())
            .collect()
    }
}
