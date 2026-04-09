import pytest
import numpy as np

from dexter import (
    LastClosedFluxSurface,
    CosHarmonic,
    Equilibrium,
    LarBfield,
    LarCurrent,
    NcBfield,
    NcCurrent,
    NcQfactor,
    NcHarmonic,
    ParabolicQfactor,
    InitialFlux,
    InitialConditions,
    IntersectParams,
    Particle,
    Perturbation,
)
from dexter.equilibrium import _POLOIDAL_TEST_NETCDF_PATH


def test_initial_flux():
    i = InitialFlux("Toroidal", 0.1)
    assert i.kind == "Toroidal"
    assert i.value == 0.1
    i = InitialFlux("Poloidal", 0.2)
    assert i.kind == "Poloidal"
    assert i.value == 0.2
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)


def test_boozer_initial_conditions():
    i = InitialConditions.boozer(
        t0=0,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    assert isinstance(i.flux0, InitialFlux)
    i = InitialConditions.boozer(
        t0=0,
        flux0=InitialFlux("Poloidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    assert i.t0 == 0
    assert isinstance(i.flux0, InitialFlux)
    assert i.theta0 == 3.14
    assert i.zeta0 == 0
    assert i.rho0 == 1e-4
    assert i.mu0 == 7e-6
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)


def test_mixed_initial_conditions():
    i = InitialConditions.mixed(
        t0=0,
        pzeta0=-0.027,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        mu0=7e-6,
    )
    assert isinstance(i.flux0, InitialFlux)
    i = InitialConditions.mixed(
        t0=0,
        pzeta0=-0.027,
        flux0=InitialFlux("Poloidal", 0.1),
        theta0=3.14,
        zeta0=0,
        mu0=7e-6,
    )
    assert i.t0 == 0
    assert i.pzeta0 == -0.027
    assert isinstance(i.flux0, InitialFlux)
    assert i.theta0 == 3.14
    assert i.zeta0 == 0
    assert i.mu0 == 7e-6
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)


def test_undefined_mixed_initial_conditions(poloidal_nc_equilibrium: Equilibrium):
    initial = InitialConditions.mixed(
        t0=0,
        pzeta0=-0.027,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        mu0=7e-6,
    )
    particle = Particle(initial)
    particle.integrate(poloidal_nc_equilibrium, (0, 1e10))
    assert particle.integration_status == "InvalidInitialConditions"


def test_intersect_params():
    IntersectParams("ConstTheta", 0, 100)
    IntersectParams("ConstZeta", 100, 10)
    with pytest.raises(TypeError):
        IntersectParams("wrong", 100, 10)  # type: ignore


def test_particle_instantiation():
    initial = InitialConditions.boozer(
        t0=0,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    particle = Particle(initial)
    assert isinstance(particle.__str__(), str)
    assert isinstance(particle.__repr__(), str)
    assert isinstance(particle.initial_conditions, InitialConditions)
    assert particle.integration_status == "Initialized"
    assert particle.initial_energy is None
    assert particle.final_energy is None
    assert particle.energy_var is None
    assert particle.steps_stored == 0
    assert particle.steps_taken == 0


# ================================================================================================


def test_mixed_toroidal_integration_default(toroidal_nc_equilibrium: Equilibrium):
    toroidal_initial = InitialConditions.mixed(
        t0=0,
        pzeta0=-0.025,
        flux0=InitialFlux("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        mu0=7e-6,
    )
    particle = Particle(toroidal_initial)
    particle.integrate(
        toroidal_nc_equilibrium,
        (0, 5),
    )
    assert particle.orbit_type == "Undefined"
    assert particle.integration_status == "Integrated"
    _check_valid_integration(particle)


def test_mixed_poloidal_integration_default(poloidal_nc_equilibrium: Equilibrium):
    poloidal_initial = InitialConditions.mixed(
        t0=0,
        pzeta0=-0.025,
        flux0=InitialFlux("Poloidal", 0.1),
        theta0=3.14,
        zeta0=0,
        mu0=7e-6,
    )
    particle = Particle(poloidal_initial)
    particle.integrate(
        poloidal_nc_equilibrium,
        (0, 3e2),
    )
    assert particle.orbit_type == "Undefined"
    assert particle.integration_status == "Integrated"
    _check_valid_integration(particle)


# ================================================================================================


class TestToroidalParticle:

    @classmethod
    def setup_class(cls) -> None:
        cls.lcfs = LastClosedFluxSurface("Toroidal", 0.45)
        cls.equilibrium = Equilibrium(
            qfactor=ParabolicQfactor(1.1, 3.9, cls.lcfs),
            current=LarCurrent(),
            bfield=LarBfield(),
            perturbation=Perturbation(
                [
                    CosHarmonic(1e-3, cls.lcfs, 1, 3, 0),
                    CosHarmonic(1e-3, cls.lcfs, 2, 3, 0),
                ]
            ),
        )
        cls.equilibrium_unperturbed = Equilibrium(
            qfactor=ParabolicQfactor(1.1, 3.9, cls.lcfs),
            current=LarCurrent(),
            bfield=LarBfield(),
            perturbation=Perturbation([]),
        )
        cls.initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Toroidal", 0.1),
            theta0=3.14,
            zeta0=0,
            rho0=1e-4,
            mu0=7e-6,
        )

    def test_integration_default(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 2e3),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_fixed_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 600),
            stepping_method=("FixedStep", 0.5),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_few_steps(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 1e20),
            max_steps=1001,
        )
        assert "TimedOut" in particle.integration_status
        _check_valid_integration(particle)

    def test_const_theta_intersection(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Toroidal", 0.1),
            theta0=0.0,
            zeta0=0.0,
            rho0=1e-6,
            mu0=0,
        )
        particle = Particle(initial)
        intersect_params = IntersectParams("ConstTheta", 1.0, 5)
        particle.intersect(self.equilibrium, intersect_params)
        assert particle.steps_stored == 5
        assert particle.integration_status == "Intersected"
        assert np.all(
            np.isclose(
                np.fmod(particle.theta_array - intersect_params.angle, 2 * np.pi), 0
            )
        )
        _check_valid_intersection(particle, intersect_params)

    def test_const_zeta_intersection(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Toroidal", 0.1),
            theta0=0.0,
            zeta0=0.0,
            rho0=1e-6,
            mu0=0,
        )
        particle = Particle(initial)
        intersect_params = IntersectParams("ConstZeta", 1.0, 5)
        particle.intersect(self.equilibrium, intersect_params)
        assert particle.steps_stored == 5
        assert particle.integration_status == "Intersected"
        assert np.all(
            np.isclose(
                np.fmod(particle.zeta_array - intersect_params.angle, 2 * np.pi), 0
            )
        )
        _check_valid_intersection(particle, intersect_params)

    def test_close(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Toroidal", 0.02),
            theta0=1,
            zeta0=0,
            rho0=1e-6,
            mu0=1e-6,
        )
        particle = Particle(initial)
        particle.close(self.equilibrium_unperturbed)
        assert particle.orbit_type == "TrappedStagnated"
        assert "ClosedPeriods" in particle.integration_status
        assert particle.omega_theta is not None
        assert particle.omega_zeta is not None
        assert particle.qkinetic is not None


class TestPoloidalParticle:

    @classmethod
    def setup_class(cls) -> None:
        path = _POLOIDAL_TEST_NETCDF_PATH
        cls.equilibrium = Equilibrium(
            qfactor=NcQfactor(path, "Steffen"),
            current=NcCurrent(path, "Steffen"),
            bfield=NcBfield(path, "Bicubic"),
            perturbation=Perturbation(
                [
                    NcHarmonic(path, "Steffen", 2, 1, "Interpolation"),
                    NcHarmonic(path, "Steffen", 3, 2, "Interpolation"),
                ]
            ),
        )
        cls.equilibrium_unperturbed = Equilibrium(
            qfactor=NcQfactor(path, "Steffen"),
            current=NcCurrent(path, "Steffen"),
            bfield=NcBfield(path, "Bicubic"),
            perturbation=Perturbation([]),
        )
        cls.initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Poloidal", 0.01),
            theta0=1.14,
            zeta0=0,
            rho0=1e-4,
            mu0=7e-6,
        )

    def test_integration_default(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 80000),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_fixed_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 300),
            stepping_method=("FixedStep", 0.1),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)
        assert particle.steps_taken == 300 / 0.1 + 2

    def test_integration_few_steps(self):
        particle = Particle(self.initial)
        particle.integrate(
            self.equilibrium,
            (0, 1e20),
            max_steps=1001,
        )
        assert "TimedOut" in particle.integration_status
        _check_valid_integration(particle)
        assert particle.steps_taken == 1001

    def test_const_theta_intersection(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Poloidal", 0.1),
            theta0=0.0,
            zeta0=0.0,
            rho0=1e-6,
            mu0=0,
        )
        particle = Particle(initial)
        intersect_params = IntersectParams("ConstTheta", 1.0, 5)
        particle.intersect(self.equilibrium, intersect_params)
        assert particle.steps_stored == 5
        assert particle.integration_status == "Intersected"
        assert np.all(
            np.isclose(
                np.fmod(particle.theta_array - intersect_params.angle, 2 * np.pi), 0
            )
        )
        _check_valid_intersection(particle, intersect_params)

    def test_const_zeta_intersection(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Poloidal", 0.1),
            theta0=0.0,
            zeta0=0.0,
            rho0=1e-6,
            mu0=0,
        )
        particle = Particle(initial)
        intersect_params = IntersectParams("ConstZeta", 1.0, 5)
        particle.intersect(self.equilibrium, intersect_params)
        assert particle.steps_stored == 5
        assert particle.integration_status == "Intersected"
        assert np.all(
            np.isclose(
                np.fmod(particle.zeta_array - intersect_params.angle, 2 * np.pi), 0
            )
        )
        _check_valid_intersection(particle, intersect_params)

    def test_close(self):
        initial = InitialConditions.boozer(
            t0=0,
            flux0=InitialFlux("Poloidal", 0.2),
            theta0=0,
            zeta0=0,
            rho0=1e-6,
            mu0=1e-6,
        )
        particle = Particle(initial)
        particle.close(self.equilibrium_unperturbed)
        assert "ClosedPeriods" in particle.integration_status
        assert particle.orbit_type == "TrappedStagnated"
        assert particle.omega_theta is not None
        assert particle.omega_zeta is not None
        assert particle.qkinetic is not None


def _check_valid_integration(particle: Particle):
    # just to keep the tests fast
    assert particle.steps_taken > 1000
    assert particle.steps_taken < 10000

    assert particle.steps_stored == particle.steps_taken

    _check_valid_time_series(particle)


def _check_valid_intersection(particle: Particle, intersect_params: IntersectParams):
    # just to keep the tests fast
    assert particle.steps_taken < 100000

    assert particle.steps_stored == intersect_params.turns


def _check_valid_time_series(particle: Particle):

    assert len(particle.t_array) == particle.steps_stored
    assert len(particle.psi_array) == particle.steps_stored
    assert len(particle.psip_array) == particle.steps_stored
    assert len(particle.theta_array) == particle.steps_stored
    assert len(particle.zeta_array) == particle.steps_stored
    assert len(particle.rho_array) == particle.steps_stored
    assert len(particle.mu_array) == particle.steps_stored
    assert len(particle.ptheta_array) == particle.steps_stored
    assert len(particle.pzeta_array) == particle.steps_stored
    assert len(particle.energy_array) == particle.steps_stored

    assert np.all(np.isfinite(particle.t_array))
    assert np.all(np.isfinite(particle.psi_array))
    assert np.all(np.isfinite(particle.psip_array))
    assert np.all(np.isfinite(particle.theta_array))
    assert np.all(np.isfinite(particle.zeta_array))
    assert np.all(np.isfinite(particle.rho_array))
    assert np.all(np.isfinite(particle.mu_array))
    assert np.all(np.isfinite(particle.ptheta_array))
    assert np.all(np.isfinite(particle.pzeta_array))
    assert np.all(np.isfinite(particle.energy_array))
