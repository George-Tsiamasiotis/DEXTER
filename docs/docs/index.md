# DEXTER - Dynamics of Experimental Toroidal Equilibrium Reconstructions

A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed equilibria from experimental data.

The data must be in [netCDF] format.

The bulk computations are implemented in *[Rust]*. The Rust code consists of standalone *[crates]* that can be used independently and *as is*.  

The Python interface directly exposes all underlying objects and routines in the form of a single python package:

```python
>>> import dexter as dx
>>>
>>> qfactor = dx.Qfactor(path="./data.nc", typ="cubic")

```

while also providing plotting methods and scripts for handling and converting netCDF files.

###  NetCDF handling

An `npz` file can be converted to a `netCDF` file with the `npz_to_netcdf` project script:

```bash
> npz_to_netcdf "data.npz" "data.nc"
```

### Plots

Plotting scripts are installed with the `dexter` python package as project scripts:

```bash
> qfactor_plot "data.nc"
> harmonic_plot "data.nc" 0 1 -t steffen
> bfield_plot --help # print help message
```

[netCDF]: https://www.unidata.ucar.edu/software/netcdf
[Rust]: https://rust-lang.org/
[crates]: https://rustup.dev/fundamentals/crates.html#crates
