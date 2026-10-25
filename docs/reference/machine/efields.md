# Efield Objects

::: dexter.EfieldObject

---

::: dexter.NoEfield
    options:
        show_bases: true
        inherited_members: false

---

::: dexter.GaussianEfield
    options:
        show_bases: true
        inherited_members: false

The GaussianEfield's profile is described by the following formulas (with $\psi$ being either the toroidal or
the poloidal flux):

$$
E(\psi) = E\cdot\exp \bigg[-\dfrac{(\psi-\mu)^2}{2\sigma^2}\bigg],
$$

$$
\Phi(\psi) = -\int_0^\psi E(\psi)d\psi = E \sigma \sqrt{\dfrac{\pi}{2}}
    \bigg[\text{erf} \bigg( \dfrac{\sqrt{2}(\psi-\mu)}{2\sigma} \bigg)
    + \text{erf} \bigg( \dfrac{\mu}{\sigma\sqrt{2}} \bigg) \bigg],
    \qquad \text{with} \quad \Phi(0) = 0
$$

$$
\dfrac{d\Phi}{d\psi} = - E(\psi),
$$

$$
\dfrac{d^2\Phi}{d\psi^2} = \dfrac{\psi - \mu}{\sigma^2}E(\psi)
$$
