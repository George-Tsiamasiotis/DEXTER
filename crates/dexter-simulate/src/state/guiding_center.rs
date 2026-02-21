use std::f64::consts::TAU;

use dexter_equilibrium::{Bfield, Current, EqError, FluxCommute, Harmonic, Perturbation, Qfactor};

use crate::Result;
use crate::particle::{Caches, EqObjects};
use crate::{InitialConditions, InitialFlux, state::FluxCoordinate};

/// State of the Guiding Center at each step
///
/// Stores all the intermediate values needed for the calculation of the final time derivatives.
///
/// Corresponds to a single specific point in configuration space, e.g. all values are calculated
/// at the same (t, ψ, θ, ζ, ρ, μ), or (t, ψp, θ, ζ, ρ, μ) point, depending on the value of the
/// `coordinate` field
#[derive(Debug, Clone)]
pub(crate) struct GCState {
    coordinate: FluxCoordinate,
    pub(crate) t: f64,

    // Dynamic variables
    pub(crate) psi: f64,
    pub(crate) psip: f64,
    pub(crate) theta: f64,
    pub(crate) zeta: f64,
    pub(crate) rho: f64,
    pub(crate) mu: f64,

    // Final time derivatives
    pub(crate) psi_dot: f64,
    pub(crate) psip_dot: f64,
    pub(crate) theta_dot: f64,
    pub(crate) zeta_dot: f64,
    pub(crate) rho_dot: f64,
    pub(crate) mu_dot: f64,

    // Modulos of angles, to avoid recalculation
    mod_theta: f64,
    mod_zeta: f64,

    pub(crate) ptheta: f64,
    pub(crate) pzeta: f64,
    pub(crate) energy: f64,

    // Equilibrium quantities
    b: f64,
    q: f64,
    g: f64,
    i: f64,
    p: f64,
    db_dflux: f64,
    db_dtheta: f64,
    db_dzeta: f64,
    dg_dflux: f64,
    di_dflux: f64,
    dp_dflux: f64,
    dp_dtheta: f64,
    dp_dzeta: f64,
    dp_dt: f64,

    /// Constants D, K, C, F, as they appear in the equations of motion.
    dterm: f64,
    kterm: f64,
    cterm: f64,
    fterm: f64,

    /// The intermediate value (μ + ρ^2Β)
    mu_par: f64,
    /// The intermediate value [dψ or dψp derivatives].
    dflux_brace: f64,
    /// The intermediate value [theta derivatives].
    dtheta_brace: f64,
    /// The intermediate value [zeta derivatives].
    dzeta_brace: f64,
    /// The intermediate value ρ*B^2/D.
    rho_bsquared_d: f64,
    /// The intermediate value g/D.
    g_over_d: f64,
    /// The intermediate value I/D.
    i_over_d: f64,
}

// Creation and evaluation
impl GCState {
    /// Creates a new `GCState` from a set of [`InitialConditions`] and evaluates it.
    pub(crate) fn new<Q, C, B, H>(
        initial: &InitialConditions,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<Self>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        let (psi, psip): (f64, f64);
        let coordinate: FluxCoordinate;
        match initial.flux0 {
            InitialFlux::Toroidal(psi0) => {
                (psi, psip) = (psi0, f64::NAN);
                coordinate = FluxCoordinate::Toroidal;
            }
            InitialFlux::Poloidal(psip0) => {
                (psi, psip) = (f64::NAN, psip0);
                coordinate = FluxCoordinate::Poloidal;
            }
        }
        Self {
            t: initial.t0,
            psi,
            psip,
            theta: initial.theta0,
            rho: initial.rho0,
            zeta: initial.zeta0,
            mu: initial.mu0,
            coordinate,
            ..Default::default()
        }
        .into_evaluated(objects, caches)
    }

    /// Performs all evaluations and calculation of intermediate quantities and final time derivatives.
    pub(crate) fn evaluate<Q, C, B, H>(
        &mut self,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        // First do all the interpolations
        self.calculate_modulos();
        self.calculate_other_flux::<Q, H>(objects.qfactor, caches)?;
        self.calculate_qfactor_quantities::<Q, H>(objects.qfactor, caches)?;
        self.calculate_current_quantities::<C, H>(objects.current, caches)?;
        self.calculate_bfield_quantities::<B, H>(objects.bfield, caches)?;
        self.calculate_perturbation_quantities::<H>(objects.perturbation, caches)?;

        // Then the intermediate quantities that only depend on the interpolated quantities
        self.calculate_canonical_momenta();
        self.calculate_capitals();
        self.calculate_mu_par();
        self.calculate_braces();
        self.calculate_extras();
        self.calculate_energy();

        // And finally the derivatives
        self.calculate_flux_dots();
        self.calculate_theta_dot();
        self.calculate_zeta_dot();
        self.calculate_rho_dot();
        self.calculate_mu_dot();
        Ok(())
    }

    /// Returns the state evaluated, consuming self.
    pub(crate) fn into_evaluated<Q, C, B, H>(
        mut self,
        objects: &EqObjects<Q, C, B, H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<Self>
    where
        Q: Qfactor + FluxCommute,
        C: Current,
        B: Bfield,
        H: Harmonic,
    {
        self.evaluate(objects, caches)?;
        Ok(self)
    }
}

/// Field calculations
impl GCState {
    fn calculate_modulos(&mut self) {
        self.mod_theta = self.theta.rem_euclid(TAU);
        self.mod_zeta = self.zeta.rem_euclid(TAU);
    }

    /// Calculates the non-coordinate flux, if it is defined.
    fn calculate_other_flux<Q, H>(
        &mut self,
        qfactor: &Q,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        H: Harmonic,
    {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.psip = match qfactor.psip_of_psi(self.psi, &mut caches.psi_acc) {
                Ok(psip) => psip,
                Err(EqError::UndefinedEvaluation(..)) => f64::NAN,
                Err(err) => return Err(err.into()),
            }
        } else {
            self.psi = match qfactor.psi_of_psip(self.psip, &mut caches.psip_acc) {
                Ok(psi) => psi,
                Err(EqError::UndefinedEvaluation(..)) => f64::NAN,
                Err(err) => return Err(err.into()),
            }
        }
        Ok(())
    }

    fn calculate_qfactor_quantities<Q, H>(
        &mut self,
        qfactor: &Q,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        Q: Qfactor + FluxCommute,
        H: Harmonic,
    {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.q = qfactor.q_of_psi(self.psi, &mut caches.psi_acc)?;
        } else {
            self.q = qfactor.q_of_psip(self.psip, &mut caches.psip_acc)?;
        };
        Ok(())
    }

    fn calculate_current_quantities<C, H>(
        &mut self,
        current: &C,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        C: Current,
        H: Harmonic,
    {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.g = current.g_of_psi(self.psi, &mut caches.psi_acc)?;
            self.i = current.i_of_psi(self.psi, &mut caches.psi_acc)?;
            self.dg_dflux = current.dg_dpsi(self.psi, &mut caches.psi_acc)?;
            self.di_dflux = current.dg_dpsi(self.psi, &mut caches.psi_acc)?;
        } else {
            self.g = current.g_of_psip(self.psip, &mut caches.psip_acc)?;
            self.i = current.i_of_psip(self.psip, &mut caches.psip_acc)?;
            self.dg_dflux = current.dg_dpsip(self.psip, &mut caches.psip_acc)?;
            self.di_dflux = current.dg_dpsip(self.psip, &mut caches.psip_acc)?;
        }
        Ok(())
    }

    #[rustfmt::skip]
    fn calculate_bfield_quantities<B, H>(
        &mut self,
        bfield: &B,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        B: Bfield,
        H: Harmonic,
    {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.b          = bfield.b_of_psi           (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dflux   = bfield.db_dpsi            (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dtheta  = bfield.db_of_psi_dtheta   (self.psi, self.mod_theta, &mut caches.psi_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dzeta   = 0.0 // Axisymmetric configuration for now
        } else {
            self.b          = bfield.b_of_psip          (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dflux   = bfield.db_dpsip           (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dtheta  = bfield.db_of_psip_dtheta  (self.psip, self.mod_theta, &mut caches.psip_acc, &mut caches.theta_acc, &mut caches.spline_cache)?;
            self.db_dzeta   = 0.0 // Axisymmetric configuration for now
        }
        Ok(())
    }

    #[rustfmt::skip]
    fn calculate_perturbation_quantities<H>(
        &mut self,
        perturbation: &Perturbation<H>,
        caches: &mut Caches<H::Cache>,
    ) -> Result<()>
    where
        H: Harmonic,
    {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.p          = perturbation.p_of_psi         (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dflux   = perturbation.dp_dpsi          (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dtheta  = perturbation.dp_of_psi_dtheta (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dzeta   = perturbation.dp_of_psi_dzeta  (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dt      = perturbation.dp_of_psi_dt     (self.psi, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
        } else {
            self.p          = perturbation.p_of_psip        (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dflux   = perturbation.dp_dpsip         (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dtheta  = perturbation.dp_of_psip_dtheta(self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dzeta   = perturbation.dp_of_psip_dzeta (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
            self.dp_dt      = perturbation.dp_of_psip_dt    (self.psip, self.mod_theta, self.mod_zeta, self.t, &mut caches.harmonic_caches)?;
        }
        Ok(())
    }

    fn calculate_canonical_momenta(&mut self) {
        self.ptheta = self.psi + self.rho * self.i;
        self.pzeta = self.rho * self.g - self.psip;
    }

    /// Calculates the matrix coefficients denoted with capital letters that appear in the
    /// perturbed equations of motion.
    ///
    /// If `ψ` is the coordinate, we must apply the chainrule in all g, I, p derivatives (with
    /// respect to `ψp`), which results in multiplying them with q.
    fn calculate_capitals(&mut self) {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.dg_dflux *= self.q;
            self.di_dflux *= self.q;
            self.dp_dflux *= self.q;
        };
        self.cterm = -1.0 + (self.rho + self.p) * self.dg_dflux + self.g * self.dp_dflux;
        self.fterm = self.q + (self.rho + self.p) * self.di_dflux + self.i * self.dp_dflux;
        self.kterm = self.g * self.dp_dtheta - self.i * self.dp_dzeta;
        self.dterm = self.g * self.fterm - self.i * self.cterm;
    }

    /// Calculates (μ + ρ^2B)
    fn calculate_mu_par(&mut self) {
        self.mu_par = self.mu + self.rho.powi(2) * self.b;
    }

    /// Calculates the brackets:
    ///     - [mu_par*dB_d<flux> + dΦ_d<flux>]
    ///     - [mu_par*dB_dtheta + dΦ_dtheta]
    /// where Φ = 0
    ///
    /// Depending on the flux coordinate, it multiplies by q-factor where needed.
    fn calculate_braces(&mut self) {
        if let FluxCoordinate::Toroidal = self.coordinate {
            self.db_dflux *= self.q;
        };
        self.dflux_brace = self.mu_par * self.db_dflux;
        self.dtheta_brace = self.mu_par * self.db_dtheta;
        self.dzeta_brace = self.mu_par * self.db_dzeta;
    }

    /// Calculates intermediate values that appear many times in the time derivatives.
    fn calculate_extras(&mut self) {
        self.rho_bsquared_d = self.rho * self.b.powi(2) / self.dterm;
        self.g_over_d = self.g / self.dterm;
        self.i_over_d = self.i / self.dterm;
    }

    fn calculate_energy(&mut self) {
        self.energy = self.energy()
    }

    fn calculate_flux_dots(&mut self) {
        let flux_dot = self.kterm * self.rho_bsquared_d - self.g_over_d * self.dtheta_brace
            + self.i_over_d * self.dzeta_brace;
        match self.coordinate {
            FluxCoordinate::Toroidal => {
                self.psi_dot = flux_dot;
                self.psip_dot = flux_dot / self.q
            }
            FluxCoordinate::Poloidal => {
                self.psi_dot = flux_dot * self.q;
                self.psip_dot = flux_dot;
            }
        }
    }

    fn calculate_theta_dot(&mut self) {
        self.theta_dot = -self.cterm * self.rho_bsquared_d + self.g_over_d * self.dflux_brace;
    }

    fn calculate_zeta_dot(&mut self) {
        self.zeta_dot = self.fterm * self.rho_bsquared_d - self.i_over_d * self.dflux_brace
    }

    fn calculate_rho_dot(&mut self) {
        self.rho_dot = self.cterm / self.dterm * self.dtheta_brace
            - self.kterm / self.dterm * self.dflux_brace
            - self.fterm / self.dterm * self.dzeta_brace
            - self.dp_dt;
    }

    fn calculate_mu_dot(&mut self) {
        self.mu_dot = 0.0;
    }

    /// Returns the Energy of the State.
    pub fn energy(&self) -> f64 {
        let parallel = self.parallel_energy();
        let perpendicular = self.perpendicular_energy();
        parallel + perpendicular
    }

    /// Returns the parallel energy of the State.
    pub fn parallel_energy(&self) -> f64 {
        // Use the ρ expression here, since the g^2 in the denominator causes numerical instability
        // on configurations with g=0.
        (self.rho * self.b).powi(2) / 2.0
    }

    /// Returns the perpendicular energy of the State.
    pub fn perpendicular_energy(&self) -> f64 {
        self.mu * self.b
    }
}

// ===============================================================================================

impl GCState {
    /// Returns a mutable reference to the flux coordinate.
    pub(crate) fn flux(&mut self) -> &mut f64 {
        match self.coordinate {
            FluxCoordinate::Toroidal => &mut self.psi,
            FluxCoordinate::Poloidal => &mut self.psip,
        }
    }

    /// Returns the array with the final time derivatives.
    pub(crate) fn dots(&self) -> [f64; 5] {
        [
            match self.coordinate {
                FluxCoordinate::Toroidal => self.psi_dot,
                FluxCoordinate::Poloidal => self.psip_dot,
            },
            self.theta_dot,
            self.zeta_dot,
            self.rho_dot,
            self.mu_dot,
        ]
    }
}

impl Default for GCState {
    /// Set all derived quantities to NaN and use the corresponding methods to set them up.
    ///
    /// If NaNs manage to propagate, it means there's something wrong.
    fn default() -> Self {
        Self {
            coordinate: Default::default(),
            t: f64::NAN,
            psi: f64::NAN,
            psip: f64::NAN,
            theta: f64::NAN,
            zeta: f64::NAN,
            rho: f64::NAN,
            mu: f64::NAN,
            psi_dot: f64::NAN,
            psip_dot: f64::NAN,
            theta_dot: f64::NAN,
            zeta_dot: f64::NAN,
            rho_dot: f64::NAN,
            mu_dot: f64::NAN,
            mod_theta: f64::NAN,
            mod_zeta: f64::NAN,
            ptheta: f64::NAN,
            pzeta: f64::NAN,
            energy: f64::NAN,
            b: f64::NAN,
            q: f64::NAN,
            g: f64::NAN,
            i: f64::NAN,
            p: f64::NAN,
            dg_dflux: f64::NAN,
            di_dflux: f64::NAN,
            db_dflux: f64::NAN,
            db_dtheta: f64::NAN,
            db_dzeta: f64::NAN,
            dp_dflux: f64::NAN,
            dp_dtheta: f64::NAN,
            dp_dzeta: f64::NAN,
            dp_dt: f64::NAN,
            dterm: f64::NAN,
            kterm: f64::NAN,
            cterm: f64::NAN,
            fterm: f64::NAN,
            mu_par: f64::NAN,
            dflux_brace: f64::NAN,
            dtheta_brace: f64::NAN,
            dzeta_brace: f64::NAN,
            rho_bsquared_d: f64::NAN,
            g_over_d: f64::NAN,
            i_over_d: f64::NAN,
        }
    }
}
