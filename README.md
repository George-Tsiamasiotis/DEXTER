# DEXTER - Dynamics of Experimental Toroidal Equilibrium Reconstructions

A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed equilibria from experimental data.

The data must be in [netCDF](https://www.unidata.ucar.edu/software/netcdf) format.

Analytical equilibria are on the (very long) [TODO](todo) list.

The bulk computations are implemented in *[Rust](https://rust-lang.org/)*. The Rust code consists of standalone *[crates](https://rustup.dev/fundamentals/crates.html#crates)* that can be used independently and *as is*.  

The Python interface directly exposes all underlying objects and routines in the form of a single python package:

``` python
>>> import dexter as dx
>>>
>>> qfactor = dx.Qfactor(path="./data.nc", typ="cubic")
```

while also providing plotting methods and scripts for handling and converting netCDF files.

###  NetCDF handling

A stub `netCDF` file can be created with the `stub_npz` project script:

``` sh
$ stub_npz ./data/stub_npz.npz
```

An `npz` file can be converted to a `netCDF` file with the `npz_to_netcdf` project script:

``` sh
$ npz_to_netcdf ./data/data.npz ./data.nc
```

### Plots

Plotting scripts are installed with the `dexter` python package as project scripts:

```sh
$ qfactor_plot ./data.nc
$ harmonic_plot .data.nc 0 1 -t steffen
$ bfield_plot --help # print help message
```
