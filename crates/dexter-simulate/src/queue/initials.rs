//! Queue Initial Conditions for batch Particle initialization.

use ndarray::Array1;

use crate::{InitialConditions, InitialFlux, Particle, SimulationError};

/// Sets of initial conditions for initializing a [`Queue`](crate::Queue`).
///
/// The [`toroidal_fluxes`] and [`poloidal_fluxes`] helper functions can be used to create the
/// initial fluxes array.
#[non_exhaustive]
#[derive(Clone)]
pub struct QueueInitialConditions {
    /// The initial times `t0`.
    t0: Array1<f64>,
    /// The initial fluxes.
    flux0: Array1<InitialFlux>,
    /// The initial thetas.
    theta0: Array1<f64>,
    /// The initial zetas.
    zeta0: Array1<f64>,
    /// The initial rhos.
    rho0: Array1<f64>,
    /// The initial mus.
    mu0: Array1<f64>,
}

impl QueueInitialConditions {
    /// Creates a new [`QueueInitialConditions`].
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// use InitialFlux::*;
    /// let initial_conditions = QueueInitialConditions::build(
    ///     &[0.0, 0.1],
    ///     &[Toroidal(0.15), Toroidal(0.3)],
    ///     &[1e-3, 2e-3],
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// # Ok::<_, SimulationError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the inputs do not have the same length.
    pub fn build(
        t0: &[f64],
        flux0: &[InitialFlux],
        theta0: &[f64],
        zeta0: &[f64],
        rho0: &[f64],
        mu0: &[f64],
    ) -> Result<Self, SimulationError> {
        let len = t0.len();
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && rho0.len() == len
            && mu0.len() == len)
        {
            return Err(SimulationError::QueueInitialConditionsMismatch);
        }
        Ok(Self {
            t0: Array1::from(t0.to_vec()),
            flux0: Array1::from(flux0.to_vec()),
            theta0: Array1::from(theta0.to_vec()),
            zeta0: Array1::from(zeta0.to_vec()),
            rho0: Array1::from(rho0.to_vec()),
            mu0: Array1::from(mu0.to_vec()),
        })
    }

    /// Creates a [`Vec<Particle>`] initialized from each set of the contained
    /// [`InitialConditions`].
    #[must_use]
    pub fn to_particles(&self) -> Vec<Particle> {
        let mut particles = Vec::with_capacity(self.len());
        for index in 0..self.len() {
            particles.push(self.particle_from_index(index));
        }
        particles
    }

    /// Creates a [`Particle`] by picking the [`InitialConditions`] at the `index` position of the arrays.
    #[must_use]
    pub fn particle_from_index(&self, index: usize) -> Particle {
        Particle::new(&self.initial_from_index(index))
    }

    /// Creates an [`InitialConditions`] from the `index` position of the arrays.
    #[must_use]
    pub(crate) fn initial_from_index(&self, index: usize) -> InitialConditions {
        InitialConditions {
            t0: self.t0[[index]],
            flux0: self.flux0[[index]],
            theta0: self.theta0[[index]],
            zeta0: self.zeta0[[index]],
            rho0: self.rho0[[index]],
            mu0: self.mu0[[index]],
        }
    }

    /// Returns the number of [`InitialConditions`] sets.
    #[must_use]
    pub fn len(&self) -> usize {
        self.t0.len()
    }

    /// Return `true` if there are no [`InitialConditions`] sets.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Helper function to create an array of `flux0`s, to be passed to [`QueueInitialConditions::build`].
///
/// # Example
/// ```
/// # use dexter_simulate::*;
/// # use ndarray::Array1;
/// let psi0s = toroidal_fluxes(&Array1::linspace(0.0, 0.45, 100).to_vec());
/// ```
#[must_use]
pub fn toroidal_fluxes(values: &[f64]) -> Array1<InitialFlux> {
    Array1::from_iter(values.iter().map(|value| InitialFlux::Toroidal(*value)))
}

/// Helper function to create an array of `flux0`s, to be passed to [`QueueInitialConditions::build`].
///
/// # Example
/// ```
/// # use dexter_simulate::*;
/// # use ndarray::Array1;
/// let psip0s = poloidal_fluxes(&Array1::linspace(0.0, 0.45, 100).to_vec());
/// ```
#[must_use]
pub fn poloidal_fluxes(values: &[f64]) -> Array1<InitialFlux> {
    Array1::from_iter(values.iter().map(|value| InitialFlux::Poloidal(*value)))
}

/// Getters.
impl QueueInitialConditions {
    /// Returns the initial times array.
    #[must_use]
    pub fn t_array(&self) -> Array1<f64> {
        self.t0.clone()
    }

    /// Returns the initial [`InitialFlux`] array.
    #[must_use]
    pub fn flux_array(&self) -> Array1<InitialFlux> {
        self.flux0.clone()
    }

    /// Returns the initial `θ` array.
    #[must_use]
    pub fn theta_array(&self) -> Array1<f64> {
        self.theta0.clone()
    }

    /// Returns the initial `ζ` array.
    #[must_use]
    pub fn zeta_array(&self) -> Array1<f64> {
        self.zeta0.clone()
    }

    /// Returns the initial `ρ` array.
    #[must_use]
    pub fn rho_array(&self) -> Array1<f64> {
        self.rho0.clone()
    }

    /// Returns the initial `μ` array.
    #[must_use]
    pub fn mu_array(&self) -> Array1<f64> {
        self.mu0.clone()
    }
}

impl std::fmt::Debug for QueueInitialConditions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("QueueInitialConditions")
            .field("Length", &self.len())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use InitialFlux::*;

    #[test]
    fn test_queue_initial_conditions_creation() {
        let p = QueueInitialConditions::build(
            &[0.0, 0.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[0.0, 1.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(p.len(), 2);
        assert!(!p.is_empty());
        let _ = format!("{:?}", p);

        let _ = QueueInitialConditions::build(
            &[0.0, 1.0, 2.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn test_queue_initial_conditions_data_extraction() {
        let p = QueueInitialConditions::build(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(p.t_array().len(), 2);
        assert_eq!(p.flux_array().len(), 2);
        assert_eq!(p.theta_array().len(), 2);
        assert_eq!(p.zeta_array().len(), 2);
        assert_eq!(p.rho_array().len(), 2);
        assert_eq!(p.mu_array().len(), 2);
    }

    #[test]
    fn test_queue_initial_conditions_to_particles() {
        let p = QueueInitialConditions::build(
            &[0.0, 1.0],
            &[Toroidal(0.1), Toroidal(0.2)],
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = p.to_particles();
        assert_eq!(particles.len(), p.len());
    }

    #[test]
    fn test_fluxes_queue_init_helpers() {
        let _: Array1<InitialFlux> = toroidal_fluxes(&[0.1, 0.2]);
        let _: Array1<InitialFlux> = poloidal_fluxes(&[0.1, 0.2]);
    }
}
