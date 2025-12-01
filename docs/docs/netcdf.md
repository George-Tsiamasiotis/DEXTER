# NetCDF Convention

Version 0.1.1

DEXTER reads equilibrium data from a [netCDF] file. The variables must follow the following conventions:

## Scalars

* `baxis`: The magnetic field strength on the magnetic axis $B_0$ in $[T]$.
* `raxis`: The horizontal position of the magnetic axis $R_0$ $[m]$.
* `zaxis`: The horizontal position of the magnetic axis $[m]$.
* `rgeo`: The geometrical axis (device major radius) in $[m]$.

!!! note

    For the normalizations, `baxis` and `raxis` should be used. They correspond to the center of the smallest flux surface. `rgeo` is characteristic of the device, while `raxis` depends on the configuration.

## Coordinates

* `psip_norm`: The poloidal flux coordinate $\psi_p$ in Normalized Units.
* `theta`: The Boozer theta coordinate $\theta_B$ in $[rads]$.
* `m`: The poloidal mode number $m$ (*index* coordinate).
* `n`: The toroidal mode number $n$ (*index* coordinate).

Not used in any calculations:

* `psi_norm`: The toroidal flux coordinate $\psi$ in Normalized Units.
* `r_norm`: The radial distance coordinate $r$ in Normalized Units.

## Variables

* `q`: The safety factor $q(\psi_p)$.
* `g_norm`: The toroidal plasma current $g(\psi_p)$ in Normalized Units.
* `i_norm`: The poloidal plasma current $I(\psi_p)$ in Normalized Units.
* `b_norm`: The magnetic field strength $B(\psi_p, \theta_B)$ in Normalized Units.
* `alphas_norm`: The harmonic amplitudes $\alpha_{m,n}(\psi_p)$ in Normalized Units.
* `phases`: The harmonic phases $\phi_{m,n}(\psi_p)$ in $[rads]$.

Original SI data (not used in any calculations):

* `r`: The radial distance coordinate $r(\psi_p)$ in $[m]$.
* `g`: The toroidal plasma current $g(\psi_p)$ in $[T \cdot m]$.
* `I`: The poloidal plasma current $I(\psi_p)$ in $[T \cdot m]$.
* `b`: The magnetic field strength $B(\psi_p, \theta_B)$ in $[T]$.
* `alphas`: The harmonic amplitudes $\alpha_{m,n}(\psi_p)$ in $[m]$.
* `rlab`: The lab horizontal coordinate $R(\psi_p, \theta_B)$ in $[m]$.
* `zlab`: The lab vertical coordinate $Z(\psi_p, \theta_B)$ in $[m]$.
* `psip`: The poloidal flux coordinate $\psi_p$ in $[T \cdot m^2]$.
* `psi`: The toroidal flux coordinate $\psi$ in $[T \cdot m^2]$.


[netCDF]: https://www.unidata.ucar.edu/software/netcdf
