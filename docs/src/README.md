# DEXTER - Dynamics of Experimental Toroidal Equilibrium Reconstructions

A code for performing calculations and simulations of particles and magnetic field lines, with reconstructed equilibria from experimental data.

The data must be in [netCDF](https://www.unidata.ucar.edu/software/netcdf) format.

Analytical equilibria are on the (very long) *TODO* list.

The bulk computations are implemented in *[Rust](https://rust-lang.org/)*. The Rust code consists of standalone *[crates](https://rustup.dev/fundamentals/crates.html#crates)* that can be used independently and *as is*.  

### Python Interface

The Python interface directly exposes all underlying objects and routines in the form of a single python package:

```python
>>> import dexter as dx 
>>> 
>>> qfactor = dx.Qfactor(path="./data.nc", typ="steffen")
>>>  
```

while also providing plotting methods and scripts for handling and converting netCDF files.


### Package scripts

Plotting scripts are installed with the `dexter` python package as project scripts:

```bash
$ qfactor_plot ./data.nc
$ harmonic_plot .data.nc 0 1 -t steffen
$ bfield_plot --help # print help message
```
