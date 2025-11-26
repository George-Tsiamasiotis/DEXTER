# Equilibrium

DEXTER uses [`interpolation`] for simulating a continuous analytic function over the equilibrium data, using the [`rsl-interpolation`] crate.

## Data interpretation

### 1D Interpolation

While still under consideration, the *[`Steffen`]* interpolation type seems to be the most suitable for 1D quantities, since it does not allow any local minima/maxima between data points, and guarantees a continuous 1st derivative.

A good comparison between the different interpolation types can be seen [`here`](<https://github.com/jfcaron3/gsl-steffen-devel/blob/steffen/interpolation/compare.pdf>)
### 2D Interpolation

At the moment, [`Bicubic Interpolation`] is the only viable option for 2D quantities.


## Perturbations

Perturbations are represented in the form:

\\[
\delta \vec B = \nabla\cdot\alpha\times \vec B
\\]

where

\\[
\alpha = \alpha(\psi_p, \theta, \zeta) = \sum\alpha_{m, n}(\psi_p)\cdot \cos\big(m\theta - n\zeta + \phi_{m,n}(\psi_p)\big)
\\]




[`interpolation`]: https://en.wikipedia.org/wiki/Interpolation
[`rsl-interpolation`]: https://docs.rs/rsl-interpolation/latest/rsl_interpolation/
[`Steffen`]: https://sourceware.org/legacy-ml/gsl-discuss/2014-q1/msg00035.html
[`Bicubic Interpolation`]: https://en.wikipedia.org/wiki/Bicubic_interpolation
