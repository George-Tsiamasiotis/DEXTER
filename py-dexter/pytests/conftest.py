import pytest
import dexter

from dexter import (
    _TEST_NETCDF_PATH as netcdf_path,
    NcCurrent,
)


@pytest.fixture(autouse=True)
def add_np(doctest_namespace):
    doctest_namespace["path"] = netcdf_path
    doctest_namespace["dex"] = dexter


@pytest.fixture(scope="session")
def nc_current() -> NcCurrent:
    """Creates an NcCurrent object from a netCDF file."""
    return NcCurrent(netcdf_path, "Steffen")
