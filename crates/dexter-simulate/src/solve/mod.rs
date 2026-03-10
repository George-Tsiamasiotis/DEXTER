//! Solver implementation.

mod rkf45;

pub(crate) use rkf45::Stepper;

#[derive(Debug, Clone)]
/// The method used to calculate the next optimal step.
pub enum SteppingMethod {
    /// Forces the step size to be small enough so that the Energy difference from step to step is
    /// under a certain threshold. The tolerances can be adjusted with the `energy_rel_tol` and
    /// `energy_abs_tol` fields.
    EnergyAdaptiveStep,
    /// Classic RK error estimation : Adjust the step size to minimize the local truncation error.
    ErrorAdaptiveStep,
    /// Fixed step size.
    FixedStep(f64),
}

/// Defines the parameters of the integration.
///
/// See [`SolverParams::default`] for the default values.
#[derive(Debug, Clone)]
pub struct SolverParams {
    /// The optimal step calculation method.
    pub method: SteppingMethod,
    /// The maximum amount of steps a particle can make before terminating its integration.
    pub max_steps: usize,
    /// The initial time step for the RKF45 adaptive step method. The value is empirical.
    pub first_step: f64,
    /// The safety factor of the solver. Should be less than 1.0.
    pub safety_factor: f64,
    /// The relative tolerance of the energy difference in every step.
    pub energy_rel_tol: f64,
    /// The absolute tolerance of the energy difference in every step.
    pub energy_abs_tol: f64,
    /// The relative tolerance of the local truncation error in every step.
    pub error_rel_tol: f64,
    /// The absolute tolerance of the local truncation error in every step.
    pub error_abs_tol: f64,
}

impl Default for SolverParams {
    fn default() -> Self {
        Self {
            method: SteppingMethod::EnergyAdaptiveStep,
            max_steps: 1_000_000,
            first_step: 1e-1,
            safety_factor: 0.9,
            energy_rel_tol: 1e-12,
            energy_abs_tol: 1e-14,
            error_rel_tol: 1e-12,
            error_abs_tol: 1e-14,
        }
    }
}

/// Defines the flux coordinate of the System.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FluxCoordinate {
    /// Use the toroidal flux `ψ` as the dynamic variable.
    #[default]
    Toroidal,
    /// Use the toroidal flux `ψp` as the dynamic variable.
    Poloidal,
}
