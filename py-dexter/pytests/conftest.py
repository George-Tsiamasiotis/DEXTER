import pytest
import dexter

from dexter import (
    _TEST_NETCDF_PATH as netcdf_path,
    NcCurrent,
    NcGeometry,
)


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace["path"] = netcdf_path
    doctest_namespace["dex"] = dexter


@pytest.fixture(scope="session")
def nc_current() -> NcCurrent:
    """Creates an NcCurrent object from the test netCDF file."""
    return NcCurrent(netcdf_path, "Steffen")


@pytest.fixture(scope="session")
def nc_geometry() -> NcGeometry:
    """Creates an NcGeometry object from the test netCDF file."""
    return NcGeometry(netcdf_path, "Steffen", "Bicubic")
