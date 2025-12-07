# Equilibrium

Handles the data extraction from [`NetCDF`] files, as well as the reconstruction of the equilibrium, providing interpolation methods for calculating all the relevant quantities.

This crate requires the [`netCDF-C`] library, which is available in most linux package managers.

`libnetcdf` can be statically linked with the `netcdf-static` feature, which is provided by the
[`netcdf crate`].

## NetCDF files

The netCDF file must follow a [`specific convention`](https://dexter.tsiamasiotis.gr/netcdf).

A *stub* netCDF file created by the `lar_netcdf.py` of the python package, which creates a netCDF file with a Large Aspect Ratio approximation equilibrium.


[`netCDF`]: https://www.unidata.ucar.edu/software/netcdf
[`netCDF-C`]: https://github.com/Unidata/netcdf-c
[`netcdf crate`]: https://github.com/georust/netcdf
