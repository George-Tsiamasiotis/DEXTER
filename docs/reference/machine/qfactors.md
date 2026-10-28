# q-factor Objects

::: dexter.QfactorObject

---

::: dexter.UnityQfactor
    options:
        show_bases: true
        inherited_members: false

---

::: dexter.ParabolicQfactor
    options:
        show_bases: true
        inherited_members: false

The ParabolicQfactor's profile is described by the following formulas:

$$
q(\psi_t) = q_{axis} + (q_{LCFS} - q_{axis})
    \bigg( \dfrac{\psi_t}{\psi_{t,LCFS}} \bigg)^2
$$

$$
\psi_p(\psi_t) = \dfrac{\psi_{t,LCFS}}{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}
    \arctan\bigg[ \dfrac{\psi_t\sqrt{q_{LCFS} - q_{axis}}}{\psi_{t,LCFS}\sqrt{q_{axis}}} \bigg]
$$

$$
\psi_t(\psi_p) = \dfrac{\psi_{t,LCFS}\sqrt{q_{axis}}}{\sqrt{q_{LCFS} - q_{axis}}}
    \tan\bigg[ \dfrac{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}{\psi_{t,LCFS}}\psi_p \bigg]
$$

$$
\dfrac{d\psi_p(\psi_t)}{d\psi_t} = ... = \dfrac{1}{q(\psi_t)} = \iota(\psi_t)
$$

$$
\dfrac{d\psi_t(\psi_p)}{d\psi_p} =
    \dfrac{q_{axis}}{\cos^2 \bigg[
    \dfrac{\sqrt{q_{axis}(q_{LCFS} - q_{axis})}}{\psi_{t,LCFS}}\psi_p
    \bigg]}
    \overset{*}{=}
    q(\psi_p)
$$

$$
q(\psi_p) = q_{axis} + q_{axis} \tan^2
    \bigg[ \dfrac{\sqrt{q_{axis}(q_{LCFS}-q_{axis})}}{\psi_{t,LCFS}} \psi_p \bigg]
$$

$$
\dfrac{dq(\psi_t)}{d\psi_t} = 2(q_{last} - q_{axis})\dfrac{\psi_t}{\psi_{t,LCFS}}
$$

$$
\dfrac{dq(\psi_p)}{d\psi_p} =
    \dfrac{2 q_{axis} \sqrt{q_{axis}(q_{last} - q_{axis})}
    \sin(\beta)}
    {\psi_{t,LCFS}\cos^3(\beta)},
    \qquad \text{where} \qquad
    \beta= \dfrac{\sqrt{q_{axis} (q_{last} - q_{axis})}\psi_p}{\psi_{t,LCFS}}
$$

$$
\iota = \dfrac{1}{q}
$$

$$
\iota' = -\dfrac{q'}{q^2}
$$

$^*$ Identity: $\dfrac{1}{\cos^2\theta} = 1 + \tan^2\theta$

---

::: dexter.NcQfactor
    options:
        show_bases: true
        inherited_members: false
