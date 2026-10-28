//! Queue Initial Conditions for batch Particle initialization.

use dexter_machine::{MagneticFluxArray, MagneticFluxKind};
use ndarray::{Array1, ArrayView1};

use crate::{CoordinateSet, InitialConditions, Particle, SimulationError};

/// Sets of initial conditions for initializing a [`Queue`](crate::Queue`).
#[non_exhaustive]
#[derive(Clone)]
pub struct QueueInitialConditions {
    /// The initial times `t0`.
    t0: Array1<f64>,
    /// The initial `fluxes`.
    flux0: MagneticFluxArray,
    /// The initial `thetas`.
    theta0: Array1<f64>,
    /// The initial `zetas`.
    zeta0: Array1<f64>,
    /// The initial `mus`.
    mu0: Array1<f64>,
    /// The initial `rhos`.
    rho0: Option<Array1<f64>>,
    /// The initial `Pzetas`.
    pzeta0: Option<Array1<f64>>,
    /// Which coordinate set were the initial conditions defined.
    coordinate_set: CoordinateSet,
}

impl QueueInitialConditions {
    /// Creates a new [`QueueInitialConditions`] in Boozer coordinates.
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the inputs do not have the same length or are empty.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// # use ndarray::*;
    /// let psit0s = MagneticFluxArray::Toroidal(array![0.15, 0.3]);
    /// let initial_conditions = QueueInitialConditions::boozer(
    ///     &[0.0, 0.1],
    ///     &psit0s,
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[1e-3, 2e-3],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// # Ok::<_, SimulationError>(())
    /// ```
    pub fn boozer(
        t0: &[f64],
        flux0: &MagneticFluxArray,
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
        if len == 0 {
            return Err(SimulationError::QueueInitialConditionsEmptyInput);
        }
        let coordinate_set = match flux0.kind() {
            MagneticFluxKind::Toroidal => CoordinateSet::BoozerToroidal,
            MagneticFluxKind::Poloidal => CoordinateSet::BoozerPoloidal,
        };

        Ok(Self {
            t0: Array1::from(t0.to_vec()),
            flux0: flux0.clone(),
            theta0: Array1::from(theta0.to_vec()),
            zeta0: Array1::from(zeta0.to_vec()),
            mu0: Array1::from(mu0.to_vec()),
            rho0: Some(Array1::from(rho0.to_vec())),
            pzeta0: None,
            coordinate_set,
        })
    }

    /// Creates a new [`QueueInitialConditions`] in Mixed coordinates.
    ///
    /// # Errors
    ///
    /// Returns a [`SimulationError`] if the inputs do not have the same length or are empty.
    ///
    /// # Example
    /// ```
    /// # use dexter_simulate::*;
    /// # use ndarray::*;
    /// let psip0s = MagneticFluxArray::Poloidal(array![0.15, 0.3]);
    /// let initial_conditions = QueueInitialConditions::mixed(
    ///     &[0.0, 0.1],
    ///     &psip0s,
    ///     &[0.0, 0.1],
    ///     &[0.0, 0.0],
    ///     &[-0.025, 0.027],
    ///     &[7e-6, 7e-6],
    /// )?;
    /// # Ok::<_, SimulationError>(())
    /// ```
    pub fn mixed(
        t0: &[f64],
        flux0: &MagneticFluxArray,
        theta0: &[f64],
        zeta0: &[f64],
        pzeta0: &[f64],
        mu0: &[f64],
    ) -> Result<Self, SimulationError> {
        let len = t0.len();
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && pzeta0.len() == len
            && mu0.len() == len)
        {
            return Err(SimulationError::QueueInitialConditionsMismatch);
        }
        if len == 0 {
            return Err(SimulationError::QueueInitialConditionsEmptyInput);
        }
        let coordinate_set = match flux0.kind() {
            MagneticFluxKind::Toroidal => CoordinateSet::MixedToroidal,
            MagneticFluxKind::Poloidal => CoordinateSet::MixedPoloidal,
        };

        Ok(Self {
            t0: Array1::from(t0.to_vec()),
            flux0: flux0.clone(),
            theta0: Array1::from(theta0.to_vec()),
            zeta0: Array1::from(zeta0.to_vec()),
            mu0: Array1::from(mu0.to_vec()),
            rho0: None,
            pzeta0: Some(Array1::from(pzeta0.to_vec())),
            coordinate_set,
        })
    }

    /// Creates a new `QueueInitialConditions` from a slice of particles by extracting their
    /// `InitialConditions`.
    ///
    /// The particles do not have to be fully initialized.
    ///
    /// # Panics
    ///
    /// Panics if *not* all the particles are defined in the same coordinate set (boozer/mixed).
    pub(crate) fn from_particles(particles: &[Particle]) -> Self {
        let mut t0 = Vec::with_capacity(particles.len());
        let mut flux0 = Vec::with_capacity(particles.len());
        let mut theta0 = Vec::with_capacity(particles.len());
        let mut zeta0 = Vec::with_capacity(particles.len());
        let mut mu0 = Vec::with_capacity(particles.len());

        // If any values are found, store them here, and the create the array
        let mut rho0 = Vec::new();
        let mut pzeta0 = Vec::new();

        for particle in particles {
            let initial = particle.initial_conditions();
            t0.push(initial.t0);
            flux0.push(initial.flux0.value());
            theta0.push(initial.theta0);
            zeta0.push(initial.zeta0);
            mu0.push(initial.mu0);
            match initial.coordinate_set() {
                CoordinateSet::BoozerToroidal | CoordinateSet::BoozerPoloidal => {
                    rho0.push(initial.rho0.expect("always some in boozer"))
                }
                CoordinateSet::MixedToroidal | CoordinateSet::MixedPoloidal => {
                    pzeta0.push(initial.pzeta0.expect("always some in mixed"))
                }
            }
        }

        // `particles` is never empty
        let flux0 = match particles[0].initial_conditions().coordinate_set() {
            CoordinateSet::BoozerToroidal | CoordinateSet::MixedToroidal => {
                MagneticFluxArray::Toroidal(Array1::from_vec(flux0))
            }
            CoordinateSet::BoozerPoloidal | CoordinateSet::MixedPoloidal => {
                MagneticFluxArray::Poloidal(Array1::from_vec(flux0))
            }
        };

        // Sanity check
        let len = t0.len();
        let other_len = rho0.len().max(pzeta0.len());
        if !(flux0.len() == len
            && theta0.len() == len
            && zeta0.len() == len
            && mu0.len() == len
            && other_len == len)
        {
            unreachable!()
        }

        Self {
            t0: Array1::from(t0),
            flux0,
            theta0: Array1::from(theta0),
            zeta0: Array1::from(zeta0),
            mu0: Array1::from(mu0),
            rho0: (!rho0.is_empty()).then(|| Array1::from(rho0)),
            pzeta0: (!pzeta0.is_empty()).then(|| Array1::from(pzeta0)),
            coordinate_set: particles[0].initial_conditions().coordinate_set(),
        }
    }

    /// Creates a [`Vec<Particle>`] initialized from each set of the contained
    /// [`InitialConditions`].
    #[must_use]
    pub fn to_particles(&self) -> Vec<Particle> {
        let mut particles = Vec::with_capacity(self.len());
        for index in 0..self.len() {
            let Some(particle) = self.particle_from_index(index) else {
                unreachable!("index is in-bounds")
            };
            particles.push(particle);
        }
        particles
    }

    /// Creates a [`Particle`] by picking the [`InitialConditions`] at the `index` position of the arrays.
    ///
    /// Returns `None` if `index` is greater than `self.len()`.
    #[must_use]
    pub fn particle_from_index(&self, index: usize) -> Option<Particle> {
        self.initial_from_index(index)
            .map(|initial| Particle::new(&initial))
    }

    /// Creates an [`InitialConditions`] from the `index` position of the arrays.
    ///
    /// Returns `None` if `index` is greater than `self.len()`.
    #[must_use]
    #[expect(clippy::missing_panics_doc, reason = "safe")]
    #[expect(clippy::unwrap_in_result, reason = "safe")]
    pub fn initial_from_index(&self, index: usize) -> Option<InitialConditions> {
        if index >= self.len() {
            return None;
        };
        use CoordinateSet::*;
        let initial = match self.coordinate_set {
            BoozerToroidal | BoozerPoloidal => InitialConditions::boozer(
                self.t0[[index]],
                self.flux0.get(index).expect("length checked"),
                self.theta0[[index]],
                self.zeta0[[index]],
                self.rho0.as_ref().expect("always some in boozer")[[index]],
                self.mu0[[index]],
            ),
            MixedToroidal | MixedPoloidal => InitialConditions::mixed(
                self.t0[[index]],
                self.flux0.get(index).expect("length checked"),
                self.theta0[[index]],
                self.zeta0[[index]],
                self.pzeta0.as_ref().expect("always some in mixed")[[index]],
                self.mu0[[index]],
            ),
        };
        Some(initial)
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

/// Getters.
impl QueueInitialConditions {
    /// Returns a view to the initial times array.
    #[must_use]
    pub fn t_array(&self) -> ArrayView1<'_, f64> {
        self.t0.view()
    }

    /// Returns a reference to the initial [`MagneticFluxArray`].
    #[must_use]
    pub fn flux_array(&self) -> &MagneticFluxArray {
        &self.flux0
    }

    /// Returns a view to the initial `θ` array.
    #[must_use]
    pub fn theta_array(&self) -> ArrayView1<'_, f64> {
        self.theta0.view()
    }

    /// Returns a view to the initial `ζ` array.
    #[must_use]
    pub fn zeta_array(&self) -> ArrayView1<'_, f64> {
        self.zeta0.view()
    }

    /// Returns a view to the initial `μ` array.
    #[must_use]
    pub fn mu_array(&self) -> ArrayView1<'_, f64> {
        self.mu0.view()
    }

    /// Returns a view to the initial `ρ` array.
    #[must_use]
    pub fn rho_array(&self) -> Option<ArrayView1<'_, f64>> {
        self.rho0.as_ref().map(|values| values.view())
    }

    /// Returns a view to the initial `Pζ` array.
    #[must_use]
    pub fn pzeta_array(&self) -> Option<ArrayView1<'_, f64>> {
        self.pzeta0.as_ref().map(|values| values.view())
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
    use crate::{IntegrationStatus, MagneticFlux};
    use ndarray::array;

    #[test]
    fn queue_boozer_initial_conditions_creation() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 0.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[0.0, 1.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(initial.len(), 2);
        assert!(!initial.is_empty());
        assert_eq!(initial.coordinate_set, CoordinateSet::BoozerToroidal);
        let _ = format!("{:?}", initial);

        let _ = QueueInitialConditions::boozer(
            &[0.0, 1.0, 2.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn queue_mixed_initial_conditions_creation() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 0.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[0.0, 1.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();
        assert_eq!(initial.len(), 2);
        assert!(!initial.is_empty());
        assert_eq!(initial.coordinate_set, CoordinateSet::MixedToroidal);
        let _ = format!("{:?}", initial);

        let _ = QueueInitialConditions::mixed(
            &[0.0, 1.0, 2.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap_err();
    }

    #[test]
    fn queue_boozer_initial_conditions_data_extraction() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 1.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(initial.t_array().len(), 2);
        assert_eq!(initial.flux_array().len(), 2);
        assert_eq!(initial.theta_array().len(), 2);
        assert_eq!(initial.zeta_array().len(), 2);
        assert_eq!(initial.mu_array().len(), 2);
        assert_eq!(initial.rho_array().unwrap().len(), 2);
        assert!(initial.pzeta_array().is_none());
    }

    #[test]
    fn queue_mixed_initial_conditions_data_extraction() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 1.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        assert_eq!(initial.t_array().len(), 2);
        assert_eq!(initial.flux_array().len(), 2);
        assert_eq!(initial.theta_array().len(), 2);
        assert_eq!(initial.zeta_array().len(), 2);
        assert_eq!(initial.mu_array().len(), 2);
        assert_eq!(initial.pzeta_array().unwrap().len(), 2);
        assert!(initial.rho_array().is_none());
    }

    #[test]
    fn queue_boozer_initial_conditions_to_particles() {
        let initial = QueueInitialConditions::boozer(
            &[0.0, 1.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = initial.to_particles();
        particles.iter().for_each(|particle| {
            assert_eq!(
                particle.integration_status(),
                IntegrationStatus::Initialized
            )
        });
        assert_eq!(particles.len(), initial.len());
    }

    #[test]
    fn queue_mixed_initial_conditions_to_particles() {
        let initial = QueueInitialConditions::mixed(
            &[0.0, 1.0],
            &MagneticFluxArray::Toroidal(array![0.1, 0.2]),
            &[1.0, 2.0],
            &[2.0, 3.0],
            &[3.0, 4.0],
            &[4.0, 5.0],
        )
        .unwrap();

        let particles: Vec<Particle> = initial.to_particles();
        particles.iter().for_each(|particle| {
            assert_eq!(
                particle.integration_status(),
                IntegrationStatus::PartiallyInitialized
            )
        });
        assert_eq!(particles.len(), initial.len());
    }

    #[test]
    fn from_boozer_particles() {
        let psit0 = MagneticFlux::Toroidal(0.001);
        let initial1 = InitialConditions::boozer(0.0, psit0, 1.0, 0.0, 0.1, 0.0);
        let initial2 = InitialConditions::boozer(0.0, psit0, 2.0, 0.0, 0.1, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }

    #[test]
    fn from_mixed_particles() {
        let psit0 = MagneticFlux::Toroidal(0.001);
        let initial1 = InitialConditions::mixed(0.0, psit0, 1.0, 0.0, 0.0, 0.0);
        let initial2 = InitialConditions::mixed(0.0, psit0, 2.0, 0.0, 0.0, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }

    #[test]
    #[should_panic]
    fn from_boozer_mixed_particles() {
        let psit0 = MagneticFlux::Toroidal(0.001);
        let initial1 = InitialConditions::mixed(0.0, psit0, 1.0, 0.0, 0.0, 0.0);
        let initial2 = InitialConditions::boozer(0.0, psit0, 2.0, 0.0, 0.1, 0.0);
        let particle1 = Particle::new(&initial1);
        let particle2 = Particle::new(&initial2);

        let initials = QueueInitialConditions::from_particles(&[particle1, particle2]);
        let queue = crate::Queue::new(&initials);

        assert_eq!(initial1.theta0, queue[0].initial_conditions().theta0);
        assert_eq!(initial2.theta0, queue[1].initial_conditions().theta0);
    }
}
