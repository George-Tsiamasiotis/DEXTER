//! Definitions of `Machine` and `MachineBuilder`.

use crate::{
    Bfield, Current, FluxCoordinateState, Geometry, MagneticFlux, MagneticFluxKind, Perturbation,
    Qfactor,
};

/// To be used when no perturbation is passed, to avoid calling `unwrap()` every time.
const ZERO_PERTURBATION: &Perturbation = &Perturbation::zero();

/// Used to create a [`Machine`].
///
/// This object only holds **references** to the machine objects.
#[derive(Debug, Clone, Copy)]
#[expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]
pub struct MachineBuilder<'obj> {
    geometry: Option<&'obj dyn Geometry>,
    qfactor: &'obj dyn Qfactor,
    current: &'obj dyn Current,
    bfield: &'obj dyn Bfield,
    perturbation: Option<&'obj Perturbation>,
}

impl<'obj> MachineBuilder<'obj> {
    /// Creates a new `MachineBuilder` from the machine's q-factor, plasma current and magnetic field.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    /// let current = LarCurrent;
    /// let bfield = LarBfield;
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn new(
        qfactor: &'obj dyn Qfactor,
        current: &'obj dyn Current,
        bfield: &'obj dyn Bfield,
    ) -> Self {
        Self {
            geometry: None,
            qfactor,
            current,
            bfield,
            perturbation: None,
        }
    }

    /// Sets the machine's [`Geometry`] object.
    ///
    /// The [`Geometry`] object is only needed for conversions to and from laboratory coordinates
    /// and is not used in any routines.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let geometry = LarGeometry::new(2.5, 1.75, 0.5);
    /// let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    /// let current = LarCurrent;
    /// let bfield = LarBfield;
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield).with_geometry(&geometry);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_geometry(mut self, geometry: &'obj dyn Geometry) -> Self {
        self.geometry = Some(geometry);
        self
    }

    /// Sets the machine's [`Perturbation`] object.
    ///
    /// If no perturbation is added, then it is set to [`Perturbation::zero()`].
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let qfactor = UnityQfactor::new(lcfs);
    /// let current = LarCurrent;
    /// let bfield = LarBfield;
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    ///
    /// let builder = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_perturbation(&perturbation);
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn with_perturbation(mut self, perturbation: &'obj Perturbation) -> Self {
        self.perturbation = Some(perturbation);
        self
    }

    /// Creates a new [`Machine`] with the Builder's configuration.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// let lcfs = MagneticFlux::Toroidal(0.45);
    /// let geometry = LarGeometry::new(2.5, 1.75, 0.5);
    /// let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    /// let current = LarCurrent;
    /// let bfield = LarBfield;
    /// let perturbation = Perturbation::new(vec![
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 2, 0.0)),
    ///     Box::new(FluteMode::new(1e-3, lcfs, 1, 3, 0.0)),
    /// ]);
    ///
    /// let machine = MachineBuilder::new(&qfactor, &current, &bfield)
    ///     .with_geometry(&geometry)
    ///     .with_perturbation(&perturbation)
    ///     .build();
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn build(self) -> Machine<'obj> {
        Machine::build(self)
    }
}

/// Container for all information about the magnetic field, geometry and perturbations of a device.
///
/// This object only holds references to the machine objects.
///
/// Should be created with a [`MachineBuilder`].
#[derive(Debug, Clone, Copy)]
#[expect(clippy::missing_docs_in_private_items, reason = "self-explanatory")]
pub struct Machine<'obj> {
    geometry: Option<&'obj dyn Geometry>,
    qfactor: &'obj dyn Qfactor,
    current: &'obj dyn Current,
    bfield: &'obj dyn Bfield,
    perturbation: &'obj Perturbation,
}

impl<'obj> Machine<'obj> {
    /// Constructs a `Machine` from a [`MachineBuilder`].
    ///
    /// # Panics
    ///
    /// This constructor panics of the `qfactor`, `current` and `bfield` objects do not have a
    /// common `good` flux coordinate.
    #[must_use]
    pub(crate) fn build(builder: MachineBuilder<'obj>) -> Self {
        let geometry = builder.geometry;
        let qfactor = builder.qfactor;
        let current = builder.current;
        let bfield = builder.bfield;
        let perturbation = builder.perturbation.unwrap_or(ZERO_PERTURBATION);

        use FluxCoordinateState::*;
        match (
            qfactor.psi_state(),
            current.psi_state(),
            bfield.psi_state(),
            qfactor.psip_state(),
            current.psip_state(),
            bfield.psip_state(),
        ) {
            (Good, Good, Good, _, _, _) | (_, _, _, Good, Good, Good) => (),
            _ => panic!(
                "'qfactor', 'current' and 'bfield' must have a common 'good' flux coordinate"
            ),
        };

        Machine {
            geometry,
            qfactor,
            current,
            bfield,
            perturbation,
        }
    }

    /// Returns the toroidal Last Closed Flux Surface.
    #[must_use]
    pub fn psi_last(&self) -> MagneticFlux {
        self.qfactor.psi_last()
    }

    /// Returns the poloidal Last Closed Flux Surface.
    #[must_use]
    pub fn psip_last(&self) -> MagneticFlux {
        self.qfactor.psip_last()
    }

    /// Returns the [`Good`] if `psi` is a good coordinate for `self.qfactor()`, `self.current()`
    /// and `self.bfield()`, else [`Bad`].
    ///
    /// [`Good`]: FluxCoordinateState::Good
    /// [`Bad`]: FluxCoordinateState::Bad
    #[must_use]
    pub fn psi_state(&self) -> FluxCoordinateState {
        use FluxCoordinateState::*;

        // `false` if even one harmonic state is `Bad`
        let mut pertubration_ok = true;
        for index in 0..self.perturbation.count() {
            if self.perturbation[index].psi_state() != Good {
                pertubration_ok = false;
                break;
            }
        }

        if (self.qfactor.psi_state() == Good)
            && (self.current.psi_state() == Good)
            && (self.bfield.psi_state() == Good)
            && pertubration_ok
        {
            Good
        } else {
            Bad
        }
    }

    /// Returns the [`Good`] if `psip` is a good coordinate for `self.qfactor()`, `self.current()`
    /// and `self.bfield()`, else [`Bad`].
    ///
    /// [`Good`]: FluxCoordinateState::Good
    /// [`Bad`]: FluxCoordinateState::Bad
    #[must_use]
    pub fn psip_state(&self) -> FluxCoordinateState {
        use FluxCoordinateState::*;

        // `false` if even one harmonic state is `Bad`
        let mut pertubration_ok = true;
        for index in 0..self.perturbation.count() {
            if self.perturbation[index].psip_state() != Good {
                pertubration_ok = false;
                break;
            }
        }

        if (self.qfactor.psip_state() == Good)
            && (self.current.psip_state() == Good)
            && (self.bfield.psip_state() == Good)
            && pertubration_ok
        {
            Good
        } else {
            Bad
        }
    }

    /// Returns [`Toroidal`] if `self.psi_state()` is [`Good`], else [`Poloidal`].
    ///
    /// Note that `psi` is checked first, so if `psi` is [`Good`], [`Toroidal`] is returned even
    /// if `psip` is [`Good`] as well.
    ///
    /// By construction of the `Machine` object, it is guaranteed that at least one of the flux
    /// will be [`Good`].
    ///
    /// [`Toroidal`]: MagneticFluxKind::Toroidal
    /// [`Poloidal`]: MagneticFluxKind::Poloidal
    /// [`Good`]: FluxCoordinateState::Good
    #[must_use]
    pub fn good_flux(&self) -> MagneticFluxKind {
        use FluxCoordinateState::*;
        match (self.psi_state(), self.psip_state()) {
            (Good, _) => MagneticFluxKind::Toroidal,
            (_, Good) => MagneticFluxKind::Poloidal,
            _ => unreachable!("states are checked at `Machine::build`"),
        }
    }

    /// Returns the the last closed flux surface of kind `kind`.
    ///
    /// Useful in routines that have already decided which flux they are going to use and need its
    /// value at the last closed flux surface.
    ///
    /// # Example
    ///
    /// ```
    /// # use dexter_machine::*;
    /// # let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
    /// # let current = LarCurrent;
    /// # let bfield = LarBfield;
    /// # let machine = MachineBuilder::new(&qfactor, &current, &bfield).build();
    /// let flux_kind = machine.good_flux();
    /// let last = machine.last(flux_kind);
    ///
    /// assert_eq!(machine.last(MagneticFluxKind::Toroidal), machine.psi_last());
    /// assert_eq!(machine.last(MagneticFluxKind::Poloidal), machine.psip_last());
    /// # Ok::<_, MachineError>(())
    /// ```
    #[must_use]
    pub fn last(&self, kind: MagneticFluxKind) -> MagneticFlux {
        match kind {
            MagneticFluxKind::Toroidal => self.psi_last(),
            MagneticFluxKind::Poloidal => self.psip_last(),
        }
    }

    /// Returns `true` if `kind` is a [`Good`][`FluxCoordinateState::Good`] magnetic flux coordinate.
    #[must_use]
    pub fn is_good_flux(&self, kind: MagneticFluxKind) -> bool {
        match kind {
            MagneticFluxKind::Toroidal => self.psi_state() == FluxCoordinateState::Good,
            MagneticFluxKind::Poloidal => self.psip_state() == FluxCoordinateState::Good,
        }
    }

    /// Returns a reference to the [`Geometry`] object.
    #[must_use]
    pub fn geometry(&self) -> Option<&dyn Geometry> {
        self.geometry
    }

    /// Returns a reference to the [`Qfactor`] object.
    #[must_use]
    pub fn qfactor(&self) -> &dyn Qfactor {
        self.qfactor
    }

    /// Returns a reference to the [`Current`] object.
    #[must_use]
    pub fn current(&self) -> &dyn Current {
        self.current
    }

    /// Returns a reference to the [`Bfield`] object.
    #[must_use]
    pub fn bfield(&self) -> &dyn Bfield {
        self.bfield
    }

    /// Returns a reference to the [`Perturbation`] object.
    #[must_use]
    pub fn perturbation(&self) -> &Perturbation {
        self.perturbation
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    #[test]
    fn machine_init() {
        let geometry = LarGeometry::new(2.5, 1.75, 0.5);
        let qfactor = UnityQfactor::new(MagneticFlux::Toroidal(0.1));
        let current = LarCurrent;
        let bfield = LarBfield;
        let perturbation = Perturbation::zero();

        let machine = MachineBuilder::new(&qfactor, &current, &bfield)
            .with_geometry(&geometry)
            .with_perturbation(&perturbation)
            .build();

        let _: Option<&dyn Geometry> = machine.geometry();
        let _: &dyn Qfactor = machine.qfactor();
        let _: &dyn Current = machine.current();
        let _: &dyn Bfield = machine.bfield();
        let _: &Perturbation = machine.perturbation();
    }
}
