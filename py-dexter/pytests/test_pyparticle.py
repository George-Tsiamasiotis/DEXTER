import numpy as np
from math import isfinite

from dexter import (
    CosHarmonic,
    InitialConditions,
    LarBfield,
    LarCurrent,
    NcBfield,
    NcCurrent,
    NcQfactor,
    NcHarmonic,
    ParabolicQfactor,
    Particle,
    perturbation,
)
from dexter._core import _PyInitialConditions
from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path


def test_initial_conditions_psi():
    i = InitialConditions(
        t0=0,
        flux0=("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    assert i.flux0 == ("Toroidal", 0.1)
    i = InitialConditions(
        t0=0,
        flux0=("Poloidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    assert isinstance(i.__str__(), str)
    assert isinstance(i.__repr__(), str)
    assert isfinite(i.t0)
    assert i.flux0 == ("Poloidal", 0.1)
    assert isfinite(i.theta0)
    assert isfinite(i.zeta0)
    assert isfinite(i.rho0)
    assert isfinite(i.mu0)


def test_particle_instantiation():
    initial = InitialConditions(
        t0=0,
        flux0=("Toroidal", 0.1),
        theta0=3.14,
        zeta0=0,
        rho0=1e-4,
        mu0=7e-6,
    )
    particle = Particle(initial)
    assert isinstance(particle.__str__(), str)
    assert isinstance(particle.__repr__(), str)
    assert isinstance(particle.initial_conditions, _PyInitialConditions)
    assert particle.integration_status == "Initialized"
    assert particle.initial_energy is None
    assert particle.final_energy is None
    assert particle.energy_var is None
    assert particle.steps_stored == 0
    assert particle.steps_taken == 0


# ================================================================================================


class TestToroidalParticle:

    @classmethod
    def setup_class(cls) -> None:
        cls.qfactor = ParabolicQfactor(1.1, 3.9, ("Toroidal", 0.45))
        cls.current = LarCurrent()
        cls.bfield = LarBfield()
        cls.per = perturbation(
            [
                CosHarmonic(1e-3, 1, 3, 0),
                CosHarmonic(1e-3, 2, 3, 0),
            ]
        )
        cls.eq = (cls.qfactor, cls.current, cls.bfield, cls.per)
        cls.initial = InitialConditions(
            t0=0,
            flux0=("Toroidal", 0.1),
            theta0=3.14,
            zeta0=0,
            rho0=1e-4,
            mu0=7e-6,
        )

    def test_integration_default(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 6),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_error_adaptive_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 600),
            method="ErrorAdaptiveStep",
            first_step=1e-5,
            error_rel_tol=1e-18,
            error_abs_tol=1e-20,
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_fixed_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 600),
            method=("FixedStep", 0.5),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)
        assert particle.steps_taken == 600 / 0.5 + 2

    def test_integration_few_steps(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 1e20),
            max_steps=1001,
        )
        assert "TimedOut" in particle.integration_status
        _check_valid_integration(particle)
        assert particle.steps_taken == 1001


class TestPoloidalParticle:

    @classmethod
    def setup_class(cls) -> None:
        cls.qfactor = NcQfactor(netcdf_path, "Steffen")
        cls.current = NcCurrent(netcdf_path, "Steffen")
        cls.bfield = NcBfield(netcdf_path, "Bicubic")
        cls.per = perturbation(
            [
                NcHarmonic(netcdf_path, "Steffen", 2, 1, "Interpolation"),
                NcHarmonic(netcdf_path, "Steffen", 3, 2, "Interpolation"),
            ]
        )
        cls.eq = (cls.qfactor, cls.current, cls.bfield, cls.per)
        cls.initial = InitialConditions(
            t0=0,
            flux0=("Poloidal", 0.1),
            theta0=3.14,
            zeta0=0,
            rho0=1e-4,
            mu0=7e-6,
        )

    def test_integration_default(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 30000),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_error_adaptive_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 6000),
            method="ErrorAdaptiveStep",
            first_step=1e-5,
            error_rel_tol=1e-18,
            error_abs_tol=1e-20,
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)

    def test_integration_fixed_step(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 300),
            method=("FixedStep", 0.1),
        )
        assert particle.integration_status == "Integrated"
        _check_valid_integration(particle)
        assert particle.steps_taken == 300 / 0.1 + 2

    def test_integration_few_steps(self):
        particle = Particle(self.initial)
        particle.integrate(
            *self.eq,
            (0, 1e20),
            max_steps=1001,
        )
        assert "TimedOut" in particle.integration_status
        _check_valid_integration(particle)
        assert particle.steps_taken == 1001


def _check_valid_integration(particle: Particle):
    # just to keep the tests fast
    assert particle.steps_taken > 1000
    assert particle.steps_taken < 10000

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
