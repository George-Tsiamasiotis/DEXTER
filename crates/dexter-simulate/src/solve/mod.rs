//! Solver implementation.

mod rkf45;

pub(crate) use rkf45::Stepper;

#[derive(Debug, Clone)]
/// The method used to calculate the next optimal step.
pub enum SteppingMethod {
    /// Forces the step size to be small enough so that the Energy difference from step to step is
    /// under a certain threshold.
    EnergyAdaptiveStep {
        /// The relative tolerance to compare the relative energy error with.
        rel_tol: f64,
        /// The absolute error tolerance. Prevents the relative energy error from becoming too
        /// small, causing the particle to get stuck.
        abs_tol: f64,
    },
    /// Classic RK error estimation: Adjust the step size to minimize the local truncation error.
    ///
    /// Note that the errors are *not* normalized, therefore the tolerances must be set according to
    /// the scale of the system's time derivatives. A good starting point is `rel_tol=1e-17` and
    /// `abs_tol=1e-19`.
    ErrorAdaptiveStep {
        /// The relative tolerance to compare the relative error with.
        rel_tol: f64,
        /// The absolute error tolerance. Prevents the relative error from becoming too
        /// small, causing the particle to get stuck.
        abs_tol: f64,
    },
    /// Fixed step size.
    FixedStep(f64),
}

impl std::fmt::Display for SteppingMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// Defines the parameters of the integration.
///
/// See [`crate::constants`] for default values.
#[derive(Debug, Clone)]
pub struct SolverParams {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method in Normalized Units. The value is
    /// empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0.
    pub safety_factor: f64,
}

impl Default for SolverParams {
    fn default() -> Self {
        #[expect(clippy::wildcard_imports, reason = "small scope")]
        use crate::constants::*;
        Self {
            method: DEFAULT_STEPPING_METHOD,
            max_steps: DEFAULT_MAX_STEPS,
            first_step: DEFAULT_FIRST_STEP,
            safety_factor: DEFAULT_SAFETY_FACTOR,
        }
    }
}

/// Defines the coordinate with respect to which the integration is performed.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IntegrationCoordinate {
    /// Use the toroidal flux `ψ` as the dynamic variable.
    #[default]
    Toroidal,
    /// Use the toroidal flux `ψp` as the dynamic variable.
    Poloidal,
}
