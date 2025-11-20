# DEXTER - Dynamics of Experimental Toroidal Equilibrium Reconstructions

A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed equilibria from experimental data.

The data must be in [netCDF](https://www.unidata.ucar.edu/software/netcdf) format.

Analytical equilibria are on the (very long) [TODO](todo) list.

The bulk computations are implemented in *[Rust](https://rust-lang.org/)*. The Rust code consists of standalone *[crates](https://rustup.dev/fundamentals/crates.html#crates)* that can be used independently and *as is*.  

The Python interface directly exposes all underlying objects and routines in the form of a single python package:

``` python
import dexter
```

while also providing plotting methods and scripts for handling and converting netCDF files.