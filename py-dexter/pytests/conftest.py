import pytest
import dexter
import matplotlib

from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path
from dexter import (
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
)


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    matplotlib.use("agg")  # Disable interactive plots
    doctest_namespace["path"] = netcdf_path
    doctest_namespace["dex"] = dexter


@pytest.fixture(scope="session")
def nc_geometry() -> NcGeometry:
    """Creates an NcGeometry object from the test netCDF file."""
    return NcGeometry(netcdf_path, "Steffen", "Bicubic")


@pytest.fixture(scope="session")
def nc_qfactor() -> NcQfactor:
    """Creates an NcQfactor object from the test netCDF file."""
    return NcQfactor(netcdf_path, "Steffen")


@pytest.fixture(scope="session")
def nc_current() -> NcCurrent:
    """Creates an NcCurrent object from the test netCDF file."""
    return NcCurrent(netcdf_path, "Steffen")


@pytest.fixture(scope="session")
def nc_bfield() -> NcBfield:
    """Creates an NcBfield object from the test netCDF file."""
    return NcBfield(netcdf_path, "Bicubic")
