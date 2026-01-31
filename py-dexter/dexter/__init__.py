from dexter._core import NcCurrent

_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/test_netcdf.nc"
_TOROIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/toroidal_test_netcdf.nc"
_POLOIDAL_TEST_NETCDF_PATH = "./crates/dexter-equilibrium/poloidal_test_netcdf.nc"

__all__ = [
    "_TOROIDAL_TEST_NETCDF_PATH",
    "_POLOIDAL_TEST_NETCDF_PATH",
    "_TEST_NETCDF_PATH",
    "NcCurrent",
]
