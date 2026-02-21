mod guiding_center;

pub(crate) use guiding_center::GCState;

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

/// Defines the flux coordinate of the System.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum FluxCoordinate {
    #[default]
    Toroidal,
    Poloidal,
}
