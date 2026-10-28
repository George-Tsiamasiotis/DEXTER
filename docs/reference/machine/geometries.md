# Geometry Objects

::: dexter.GeometryObject

---

::: dexter.LarGeometry
    options:
        show_bases: true
        inherited_members: false

A LAR configuration's geometry is described by the following formulas:

$$
r(\psi_t) = \sqrt{2\psi_t}
$$

$$
\psi_t(r) = r^2/2
$$

$$
R(\psi_t, \theta) = R_{geo} + R_{geo}\sqrt{2\psi_t}\cos\theta
$$

$$
Z(\psi_t, \theta) = R_{geo}\sqrt{2\psi_t}\sin\theta
$$

---

::: dexter.NcGeometry
    options:
        show_bases: true
        inherited_members: false
