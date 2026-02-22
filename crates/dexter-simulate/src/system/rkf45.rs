/// Runge-kutta-Fehlberg method of order 4(5).
use self::tableau::*;
use super::{StepperParams, SteppingMethod};
use crate::{
    Result, SimulationError,
    particle::{Caches, EqObjects},
    state::GCState,
};
use dexter_equilibrium::{Bfield, Current, FluxCommute, Harmonic, Qfactor};

/// Runge-kutta-Fehlberg method coefficients (Wikipedia)
mod tableau {

    // Fifth order
    pub(super) const B1: f64 = 25.0 / 216.0;
    pub(super) const B2: f64 = 0.0;
    pub(super) const B3: f64 = 1408.0 / 2565.0;
    pub(super) const B4: f64 = 2197.0 / 4104.0;
    pub(super) const B5: f64 = -1.0 / 5.0;
    pub(super) const B6: f64 = 0.0;

    // Embedded
    pub(super) const B1E: f64 = 16.0 / 135.0;
    pub(super) const B2E: f64 = 0.0;
    pub(super) const B3E: f64 = 6656.0 / 12825.0;
    pub(super) const B4E: f64 = 28561.0 / 56430.0;
    pub(super) const B5E: f64 = -9.0 / 50.0;
    pub(super) const B6E: f64 = 2.0 / 55.0;

    // pub(super) const C1: f64 = 0.0; Always equals to 0.
    pub(super) const C2: f64 = 1.0 / 4.0;
    pub(super) const C3: f64 = 3.0 / 8.0;
    pub(super) const C4: f64 = 12.0 / 13.0;
    pub(super) const C5: f64 = 1.0;
    pub(super) const C6: f64 = 1.0 / 2.0;

    pub(super) const A21: f64 = 1.0 / 4.0;

    pub(super) const A31: f64 = 3.0 / 32.0;
    pub(super) const A32: f64 = 9.0 / 32.0;

    pub(super) const A41: f64 = 1932.0 / 2197.0;
    pub(super) const A42: f64 = -7200.0 / 2197.0;
    pub(super) const A43: f64 = 7296.0 / 2197.0;

    pub(super) const A51: f64 = 439.0 / 216.0;
    pub(super) const A52: f64 = -8.0;
    pub(super) const A53: f64 = 3680.0 / 513.0;
    pub(super) const A54: f64 = -845.0 / 4104.0;

    pub(super) const A61: f64 = -8.0 / 27.0;
    pub(super) const A62: f64 = 2.0;
    pub(super) const A63: f64 = -3544.0 / 2565.0;
    pub(super) const A64: f64 = 1859.0 / 4104.0;
    pub(super) const A65: f64 = -11.0 / 40.0;
}

#[derive(Debug)]
pub(crate) struct Stepper {
    k1: [f64; 5],
    k2: [f64; 5],
    k3: [f64; 5],
    k4: [f64; 5],
    k5: [f64; 5],
    k6: [f64; 5],
    state1: GCState,
    state2: GCState,
    state3: GCState,
    state4: GCState,
    state5: GCState,
    state6: GCState,
    weights: [f64; 5],
    pub(crate) errors: [f64; 5],
}

impl Stepper {
    /// Initializes a [`Stepper`] from an *evaluated* state
    pub(crate) fn new(state: &GCState) -> Self {
        Self {
            state1: state.clone(),
            ..Default::default()
        }
    }

    /// Calculates all intermediate [`SystemState`]s, coefficients, weights and errors.
    pub(crate) fn start<C, Q, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        self.calculate_k1();
        self.calculate_state_k2(h, objects, caches)?;
        self.calculate_state_k3(h, objects, caches)?;
        self.calculate_state_k4(h, objects, caches)?;
        self.calculate_state_k5(h, objects, caches)?;
        self.calculate_state_k6(h, objects, caches)?;
        self.calculate_embedded_weights();
        self.calculate_errors();
        // Check if any NaNs where produced/propagated during the evaluations. `self.weights`
        // depends on all previous computations, so if any NaN appears, it will show up here.
        // If we dont check it here, next_optimal_step() will return NaN, which would be much
        // harder to debug.
        self.weights
            .iter()
            .try_for_each(|val| match val.is_finite() {
                false => {
                    // This would be a truly peculiar situation, and might be missed when dealing
                    // with a lot of particles, so lets print it as well. Still unsure if this
                    // should be a panic or simply discard the particle.
                    dbg!(&self);
                    eprintln!("{:?}", SimulationError::StepperNan);
                    Err(SimulationError::StepperNan)
                }
                true => Ok(()),
            })
    }

    pub(crate) fn calculate_k1(&mut self) {
        self.k1 = self.state1.dots()
    }

    pub(crate) fn calculate_state_k2<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let coef = [
            A21 * self.k1[0],
            A21 * self.k1[1],
            A21 * self.k1[2],
            A21 * self.k1[3],
            A21 * self.k1[4],
        ];

        self.state2.coordinate = self.state1.coordinate;
        self.state2.t = self.state1.t + C2 * h;
        *self.state2.flux() = *self.state1.flux() + coef[0] * h;
        self.state2.theta = self.state1.theta + coef[1] * h;
        self.state2.zeta = self.state1.zeta + coef[2] * h;
        self.state2.rho = self.state1.rho + coef[3] * h;
        self.state2.mu = self.state1.mu + coef[4] * h;
        self.state2.evaluate(objects, caches)?;
        self.k2 = self.state2.dots();
        Ok(())
    }

    pub(crate) fn calculate_state_k3<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let coef = [
            A31 * self.k1[0] + A32 * self.k2[0],
            A31 * self.k1[1] + A32 * self.k2[1],
            A31 * self.k1[2] + A32 * self.k2[2],
            A31 * self.k1[3] + A32 * self.k2[3],
            A31 * self.k1[4] + A32 * self.k2[4],
        ];

        self.state3.coordinate = self.state1.coordinate;
        self.state3.t = self.state1.t + C3 * h;
        *self.state3.flux() = *self.state1.flux() + coef[0] * h;
        self.state3.theta = self.state1.theta + coef[1] * h;
        self.state3.zeta = self.state1.zeta + coef[2] * h;
        self.state3.rho = self.state1.rho + coef[3] * h;
        self.state3.mu = self.state1.mu + coef[4] * h;
        self.state3.evaluate(objects, caches)?;
        self.k3 = self.state3.dots();
        Ok(())
    }

    pub(crate) fn calculate_state_k4<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let coef = [
            A41 * self.k1[0] + A42 * self.k2[0] + A43 * self.k3[0],
            A41 * self.k1[1] + A42 * self.k2[1] + A43 * self.k3[1],
            A41 * self.k1[2] + A42 * self.k2[2] + A43 * self.k3[2],
            A41 * self.k1[3] + A42 * self.k2[3] + A43 * self.k3[3],
            A41 * self.k1[4] + A42 * self.k2[4] + A43 * self.k3[4],
        ];

        self.state4.coordinate = self.state1.coordinate;
        self.state4.t = self.state1.t + C4 * h;
        *self.state4.flux() = *self.state1.flux() + coef[0] * h;
        self.state4.theta = self.state1.theta + coef[1] * h;
        self.state4.zeta = self.state1.zeta + coef[2] * h;
        self.state4.rho = self.state1.rho + coef[3] * h;
        self.state4.mu = self.state1.mu + coef[4] * h;
        self.state4.evaluate(objects, caches)?;
        self.k4 = self.state4.dots();
        Ok(())
    }

    pub(crate) fn calculate_state_k5<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let coef = [
            A51 * self.k1[0] + A52 * self.k2[0] + A53 * self.k3[0] + A54 * self.k4[0],
            A51 * self.k1[1] + A52 * self.k2[1] + A53 * self.k3[1] + A54 * self.k4[1],
            A51 * self.k1[2] + A52 * self.k2[2] + A53 * self.k3[2] + A54 * self.k4[2],
            A51 * self.k1[3] + A52 * self.k2[3] + A53 * self.k3[3] + A54 * self.k4[3],
            A51 * self.k1[4] + A52 * self.k2[4] + A53 * self.k3[4] + A54 * self.k4[4],
        ];

        self.state5.coordinate = self.state1.coordinate;
        self.state5.t = self.state1.t + C5 * h;
        *self.state5.flux() = *self.state1.flux() + coef[0] * h;
        self.state5.theta = self.state1.theta + coef[1] * h;
        self.state5.zeta = self.state1.zeta + coef[2] * h;
        self.state5.rho = self.state1.rho + coef[3] * h;
        self.state5.mu = self.state1.mu + coef[4] * h;
        self.state5.evaluate(objects, caches)?;
        self.k5 = self.state5.dots();
        Ok(())
    }

    pub(crate) fn calculate_state_k6<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        #[rustfmt::skip]
            let coef = [
                A61 * self.k1[0] + A62 * self.k2[0] + A63 * self.k3[0] + A64 * self.k4[0] + A65 * self.k5[0],
                A61 * self.k1[1] + A62 * self.k2[1] + A63 * self.k3[1] + A64 * self.k4[1] + A65 * self.k5[1],
                A61 * self.k1[2] + A62 * self.k2[2] + A63 * self.k3[2] + A64 * self.k4[2] + A65 * self.k5[2],
                A61 * self.k1[3] + A62 * self.k2[3] + A63 * self.k3[3] + A64 * self.k4[3] + A65 * self.k5[3],
                A61 * self.k1[4] + A62 * self.k2[4] + A63 * self.k3[4] + A64 * self.k4[4] + A65 * self.k5[4],
            ];

        self.state6.coordinate = self.state1.coordinate;
        self.state6.t = self.state1.t + C6 * h;
        *self.state6.flux() = *self.state1.flux() + coef[0] * h;
        self.state6.theta = self.state1.theta + coef[1] * h;
        self.state6.zeta = self.state1.zeta + coef[2] * h;
        self.state6.rho = self.state1.rho + coef[3] * h;
        self.state6.mu = self.state1.mu + coef[4] * h;
        self.state6.evaluate(objects, caches)?;
        self.k6 = self.state6.dots();
        Ok(())
    }

    pub(crate) fn calculate_embedded_weights(&mut self) {
        for i in 0..self.weights.len() {
            self.weights[i] = B1 * self.k1[i]
                + B2E * self.k2[i]
                + B3E * self.k3[i]
                + B4E * self.k4[i]
                + B5E * self.k5[i]
                + B6E * self.k6[i];
        }
    }
    pub(crate) fn calculate_errors(&mut self) {
        for i in 0..self.errors.len() {
            self.errors[i] = (B1 - B1E) * self.k1[i]
                + (B2 - B2E) * self.k2[i]
                + (B3 - B3E) * self.k3[i]
                + (B4 - B4E) * self.k4[i]
                + (B5 - B5E) * self.k5[i]
                + (B6 - B6E) * self.k6[i];
        }
    }

    pub(crate) fn calculate_optimal_step(&self, h: f64, params: &impl StepperParams) -> f64 {
        match params.method() {
            SteppingMethod::EnergyAdaptiveStep => self.energy_method_optimal_step(h, params),
            SteppingMethod::ErrorAdaptiveStep => self.error_method_optimal_step(h, params),
            SteppingMethod::FixedStep(stepsize) => *stepsize,
        }
    }

    /// Adjust the error by calculating the relative difference in the energy at every step.
    ///
    /// Source:
    /// https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ss2017/numerische_Methoden_fuer_komplexe_Systeme_II/rkm-1.pdf
    fn energy_method_optimal_step(&self, h: f64, config: &impl StepperParams) -> f64 {
        let initial_energy = self.state1.energy();
        let final_energy = self.state6.energy();
        // When the energy diff happens to be smaller than REL_TOL, the optimal step keeps getting
        // smaller due to the `REL_TOL/energy_diff` factor, so we need to bound it
        let energy_diff = ((initial_energy - final_energy) / initial_energy)
            .abs()
            .max(config.energy_abs_tol());
        let exp = if energy_diff >= config.energy_rel_tol() {
            0.2
        } else {
            0.25
        };
        config.safety_factor() * h * (config.energy_rel_tol() / energy_diff).powf(exp)
    }

    /// Adjust the error by calculating the local truncation error.
    ///
    /// Source:
    /// https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ss2017/numerische_Methoden_fuer_komplexe_Systeme_II/rkm-1.pdf
    fn error_method_optimal_step(&self, h: f64, config: &impl StepperParams) -> f64 {
        // Using the max error vs each variable's error is equivalent.

        // The only way this could fail was if `self.errors` contained any non-finite values,
        // which is already checked at the end of `start()`.
        let mut max_error = self
            .errors
            .iter()
            .max_by(|a, b| a.abs().total_cmp(&b.abs()))
            .copied()
            .expect("only finite values here");

        // When all errors happen to be smaller than REL_TOL, the optimal step keeps getting
        // smaller due to the `REL_TOL/max_error` factor, so we need to bound it
        max_error = max_error.max(config.error_abs_tol());

        // 0.2 = 1/(p+1), where p the order
        let exp = if max_error >= config.error_rel_tol() {
            0.2
        } else {
            0.25
        };
        config.safety_factor() * h * (config.error_rel_tol() / max_error).powf(exp)
    }

    pub(crate) fn next_state<Q, C, B, H>(
        &mut self,
        h: f64,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<GCState>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        {
            let mut next = GCState::default();
            next.coordinate = self.state1.coordinate;
            next.t = self.state1.t + h;
            *next.flux() = *self.state1.flux() + h * self.weights[0];
            next.theta = self.state1.theta + h * self.weights[1];
            next.zeta = self.state1.zeta + h * self.weights[2];
            next.rho = self.state1.rho + h * self.weights[3];
            next.mu = self.state1.mu + h * self.weights[4];
            next.into_evaluated(objects, caches)
        }
    }
}

impl Default for Stepper {
    fn default() -> Self {
        Self {
            k1: [f64::NAN; 5],
            k2: [f64::NAN; 5],
            k3: [f64::NAN; 5],
            k4: [f64::NAN; 5],
            k5: [f64::NAN; 5],
            k6: [f64::NAN; 5],
            weights: [f64::NAN; 5],
            errors: [f64::NAN; 5],
            state1: GCState::default(),
            state2: GCState::default(),
            state3: GCState::default(),
            state4: GCState::default(),
            state5: GCState::default(),
            state6: GCState::default(),
        }
    }
}
