# Project Scripts

`DEXTER` comes with some ready-to-use project scripts. These scripts become available as commands when the
`dexter` package is installed through pip. Running `<command> -h` prints info on how to use them.

- `bfield_plot`: Creates contour plots of the magnetic field and its derivatives from a `netCDF` file.

# Tools

- `./tools/create_test_netcdf.py`: Creates the 3 different `netCDF` files that are needed for unit testing.

- `./tools/create_test_netcdf.sh`: Creates the 3 test `netCDF` files and places them at the correct locations.

- `./tools/falcon_to_dexter.py`: Converts a `falcon` `netCDF` to a `DEXTER` one. Only changes the variable names.
