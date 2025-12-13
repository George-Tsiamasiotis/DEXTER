import pytest

import dexter
from dexter import (
    _LAR_NETCDF_PATH as netcdf_path,
    NcGeometry,
    NcQfactor,
    UnityQfactor,
    NcCurrent,
    LarCurrent,
    NcBfield,
    NcHarmonic,
    NcPerturbation,
)


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace["path"] = netcdf_path
    doctest_namespace["dex"] = dexter


@pytest.fixture(scope="session")
def nc_geometry() -> NcGeometry:
    """Creates an NcGeometry object from a netCDF file."""
    return NcGeometry(netcdf_path, "Steffen", "Bicubic")


@pytest.fixture(scope="session")
def nc_qfactor() -> NcQfactor:
    """Creates a Qfactor object from a netCDF file."""
    return NcQfactor(netcdf_path, "Steffen")


@pytest.fixture(scope="session")
def unity_qfactor() -> UnityQfactor:
    """Creates a Unity Qfactor object."""
    return UnityQfactor()


@pytest.fixture(scope="session")
def nc_current() -> NcCurrent:
    """Creates an NcCurrent object from a netCDF file."""
    return NcCurrent(netcdf_path, "Steffen")


@pytest.fixture(scope="session")
def lar_current() -> LarCurrent:
    """Creates an LarCurrent object."""
    return LarCurrent()


@pytest.fixture(scope="session")
def nc_bfield() -> NcBfield:
    """Creates a NcBfield object from a netCDF file."""
    return NcBfield(netcdf_path, "Bicubic")


@pytest.fixture(scope="session")
def nc_harmonic1() -> NcHarmonic:
    """Creates a NcHarmonic object from a netCDF file."""
    return NcHarmonic(netcdf_path, "Akima", m=2, n=1)


@pytest.fixture(scope="session")
def nc_harmonic2() -> NcHarmonic:
    """Creates a NcHarmonic object from a netCDF file."""
    return NcHarmonic(netcdf_path, "Akima", m=3, n=2)


@pytest.fixture(scope="session")
def nc_perturbation(
    nc_harmonic1: NcHarmonic, nc_harmonic2: NcHarmonic
) -> NcPerturbation:
    """Creates a Perturbation object with the 2 fixture harmonics."""
    return NcPerturbation(harmonics=[nc_harmonic1, nc_harmonic2])
