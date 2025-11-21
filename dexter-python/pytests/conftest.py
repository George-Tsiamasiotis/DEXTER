import pytest

import dexter as dx


netcdf_path = "./data.nc"


@pytest.fixture(scope="session")
def qfactor():
    """Creates a Qfactor object from a netCDF file."""
    return dx.Qfactor(netcdf_path, "akima")


@pytest.fixture(scope="session")
def currents():
    """Creates a Current object from a netCDF file."""
    return dx.Currents(netcdf_path, "akima")


@pytest.fixture(scope="session")
def bfield():
    """Creates a Bfield object from a netCDF file."""
    return dx.Bfield(netcdf_path, "bicubic")


@pytest.fixture(scope="session")
def harmonic1():
    """Creates a Harmonic object from a netCDF file."""
    return dx.Harmonic(netcdf_path, "akima", m=1, n=2)


@pytest.fixture(scope="session")
def harmonic2():
    """Creates a Harmonic object from a netCDF file."""
    return dx.Harmonic(netcdf_path, "akima", m=1, n=3)


@pytest.fixture(scope="session")
def perturbation(harmonic1, harmonic2):
    """Creates a Perturbation object with the 2 fixture harmonics."""
    return dx.Perturbation(harmonics=[harmonic1, harmonic2])
