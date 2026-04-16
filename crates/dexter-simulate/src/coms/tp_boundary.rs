//! Representation of the Trapped-Passing boundary curves on the `(E, Pζ, μ=const)` space.

use std::f64::consts::PI;

use dexter_equilibrium::{Bfield, FluxCommute, FluxCoordinateState, Qfactor};
use ndarray::Array1;
use rsl_interpolation::{Accelerator, Cache};

use crate::constants::TRAPPED_PASSING_BOUNDARY_DENSITY;

/// Representation of the Trapped-Passing boundary curves on the `(E, Pζ, μ=const)` space.
///
/// The Trapped-Passing boundary is defined by the curves:
///
/// 1. E = μB(ψp, 0)
/// 2. E = μB(ψp, π)
///
/// where `Pζ = -ψp`. Therefore, we can rewrite them in the `E-Pζ` plane.
///
/// 1. E(Pζ) = μB(-Pζ, 0)
/// 2. E(Pζ) = μB(-Pζ, π)
///
/// The peaks of the 2 wall parabolas and the axis parabola are always at the `Pζ/ψp_last = -1` and
/// `Pζ/ψp_last = 0` respectively. Therefore, the two curves lie in the `[-ψp_last, 0]` interval.
///
/// # Note
///
/// If the equilibrium has a `good` ψp coordinate, then it is used for all calculations, since
/// its faster. This also results in an equispaced [`TrappedPassingBoundary::pzeta_interval`].
///
/// If ψ is the only `good` coordinate, some extra conversions are necessary and the
/// `pzeta_interval` is no longer a linspace.
///
#[derive(Debug, Clone)]
pub struct TrappedPassingBoundary {
    /// The `Pζ = [-ψp_last, 0]` interval array.
    pub(crate) pzeta_interval: Array1<f64>,
    /// The curve corresponding to the lower part of the boundary, defined by `θ=0` (eq. 1).
    pub(crate) lower: Array1<f64>,
    /// The curve corresponding to the upper part of the boundary, defined by `θ=π` (eq. 2).
    pub(crate) upper: Array1<f64>,
}

impl TrappedPassingBoundary {
    /// Calculates the two curves.
    pub fn new<Q, B>(qfactor: &Q, bfield: &B, mu: f64) -> Self
    where
        Q: Qfactor + FluxCommute,
        B: Bfield,
    {
        // try `psip` first, since its faster.
        if bfield.psip_state() == FluxCoordinateState::Good {
            Self::from_good_psip(qfactor, bfield, mu)
        } else if bfield.psi_state() == FluxCoordinateState::Good {
            Self::from_good_psi(qfactor, bfield, mu)
        } else {
            unreachable!()
        }
    }

    /// Returns the [`TrappedPassingBoundary`]'s `Pζ = [-ψp_last, 0]` interval array.
    #[must_use]
    pub fn pzeta_interval(&self) -> Array1<f64> {
        self.pzeta_interval.clone()
    }

    /// Returns the [`TrappedPassingBoundary`]'s upper curve.
    #[must_use]
    pub fn upper(&self) -> Array1<f64> {
        self.upper.clone()
    }

    /// Returns the [`TrappedPassingBoundary`]'s lower curve.
    #[must_use]
    pub fn lower(&self) -> Array1<f64> {
        self.lower.clone()
    }
}

impl TrappedPassingBoundary {
    /// Creates the boundary in a `Good ψp` equilibrium.
    ///
    /// This method cannot fail since `FluxState` is already checked and `ψp` is always in-bounds.
    fn from_good_psip<Q, B>(qfactor: &Q, bfield: &B, mu: f64) -> Self
    where
        Q: Qfactor + FluxCommute,
        B: Bfield,
    {
        let psip_last = qfactor.psip_last();
        let psip_interval = Array1::linspace(0.0, psip_last, TRAPPED_PASSING_BOUNDARY_DENSITY);

        let mut psip_acc = Accelerator::new();
        let mut theta_acc = Accelerator::new();
        let mut cache = Cache::new();

        let lower = mu
            * psip_interval.mapv(|psip| {
                bfield
                    .b_of_psip(psip, 0.0, &mut psip_acc, &mut theta_acc, &mut cache)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });
        let upper = mu
            * psip_interval.mapv(|psip| {
                bfield
                    .b_of_psip(psip, PI, &mut psip_acc, &mut theta_acc, &mut cache)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });
        Self {
            pzeta_interval: -&psip_interval,
            upper,
            lower,
        }
    }

    /// Creates the boundary in a `Good ψ` equilibrium.
    ///
    /// This method cannot fail since `FluxState` is already checked and `ψp/ψ` is always in-bounds.
    ///
    /// NOTE:
    /// We define the `psi_interval` first and the `psip_interval` second, since it is not
    /// guaranteed that `qfactor` defines `ψ(ψp)`. This results in a non-linspace `pzeta_interval`,
    /// but the values are correct.
    fn from_good_psi<Q, B>(qfactor: &Q, bfield: &B, mu: f64) -> Self
    where
        B: Bfield,
        Q: Qfactor + FluxCommute,
    {
        let psi_last = qfactor.psi_last();
        let psi_interval = Array1::linspace(0.0, psi_last, TRAPPED_PASSING_BOUNDARY_DENSITY);

        let mut psi_acc = Accelerator::new();
        let mut theta_acc = Accelerator::new();
        let mut cache = Cache::new();

        let lower = mu
            * psi_interval.mapv(|psi| {
                bfield
                    .b_of_psi(psi, 0.0, &mut psi_acc, &mut theta_acc, &mut cache)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });
        let upper = mu
            * psi_interval.mapv(|psi| {
                bfield
                    .b_of_psi(psi, PI, &mut psi_acc, &mut theta_acc, &mut cache)
                    .expect("-pzeta=psip is always inbound and evaluation is defined")
            });

        let psip_interval = psi_interval.mapv(|psi| {
            qfactor
                .psip_of_psi(psi, &mut psi_acc)
                .expect("psi is always inbound and evaluation is defined")
        });

        Self {
            pzeta_interval: -&psip_interval,
            upper,
            lower,
        }
    }
}
