import pytest

from dexter import Qfactor, Currents, Bfield, Harmonic, Perturbation
from dexter import InitialConditions


netcdf_path = "./data.nc"


@pytest.fixture(scope="session")
def qfactor() -> Qfactor:
    """Creates a Qfactor object from a netCDF file."""
    return Qfactor(netcdf_path, "akima")


@pytest.fixture(scope="session")
def currents() -> Currents:
    """Creates a Current object from a netCDF file."""
    return Currents(netcdf_path, "akima")


@pytest.fixture(scope="session")
def bfield() -> Bfield:
    """Creates a Bfield object from a netCDF file."""
    return Bfield(netcdf_path, "bicubic")


@pytest.fixture(scope="session")
def harmonic1() -> Harmonic:
    """Creates a Harmonic object from a netCDF file."""
    return Harmonic(netcdf_path, "akima", m=1, n=2)


@pytest.fixture(scope="session")
def harmonic2() -> Harmonic:
    """Creates a Harmonic object from a netCDF file."""
    return Harmonic(netcdf_path, "akima", m=1, n=3)


@pytest.fixture(scope="session")
def perturbation(harmonic1: Harmonic, harmonic2: Harmonic) -> Perturbation:
    """Creates a Perturbation object with the 2 fixture harmonics."""
    return Perturbation(harmonics=[harmonic1, harmonic2])


# ==========================================================================


@pytest.fixture(scope="session")
def initial_conditions(qfactor: Qfactor) -> InitialConditions:
    """Creates an InitialConditions set."""
    return InitialConditions(
        time0=0,
        theta0=0,
        psip0=qfactor.psip_wall,
        rho0=1e-3,
        zeta0=0,
        mu=0,
    )
