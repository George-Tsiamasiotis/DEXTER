"""Defines the different efield objects as wrappers over `_PyEfield`."""

from typing import TypeAlias

from dexter._core import _PyEfield

from dexter.machine.base import Efield, MachineObject
from dexter.machine.flux import MagneticFlux
from pint.facets.plain import PlainQuantity


class NoEfield(MachineObject, Efield):
    r"""No electric field.

    You shouldn't have to use this type; It is created when no `efield` parameter is passed to
    [`Machine`][dexter.Machine].

    No $\psi/\psi_p$ bounds checks are performed in evaluations.

    Example
    -------
    ``` py
    >>> efield = dex.NoEfield()

    ```
    """

    _r: _PyEfield

    def __init__(self) -> None:
        self._r = _PyEfield.build_no_efield()
        MachineObject.__init__(self)
        Efield.__init__(self)


class GaussianEfield(MachineObject, Efield):
    r"""Radial-Gaussian electric field.

    Parameters
    ----------
    strength
        The electric field strength at the curve's peak. If `float` then it must be in *Normalized
        Units*. If [`Quantity`](https://pint.readthedocs.io/en/stable/user/defining-quantities.html)
        then it must have the dimensionality of $[V/m]$. Defining it through an existing
        [`Machine's`][dexter.Machine] [`quantity`][dexter.Machine.quantity] method is necessary to properly
        convert to Normalized Units (see example). The `machine`'s [`GeometryObject`][dexter.GeometryObject]
        *must* be defined.
    mean
        The magnetic flux surface where the peak is located.
    width
        The *standard deviation* or *root mean square width* of the 'bell'.

    Attributes
    ----------
    strength
        The electric field strength at the curve's peak, **in Normalized Units**.
    mean
        The magnetic flux surface where the peak is located.
    width
        The *standard deviation* or *root mean square width* of the 'bell'.

    Example
    -------
    ``` py title="Electric Field in Normalized Units"
    >>> lcfs = dex.MagneticFlux.Toroidal(0.05)
    >>>
    >>> strength = 1e-4
    >>> mean = lcfs / 2
    >>> width = lcfs.value / 50
    >>>
    >>> efield = dex.GaussianEfield(strength, mean, width)

    ```

    ``` py title="Electric Field from `machine`"
    >>> geometry=dex.LarGeometry(baxis=2, raxis=2.75, rlast=0.5)
    >>> lcfs = geometry.psi_last
    >>> machine = dex.Machine(
    ...     geometry=geometry,
    ...     qfactor=dex.UnityQfactor(lcfs),
    ...     current=dex.LarCurrent(),
    ...     bfield=dex.LarBfield(),
    ...     perturbation=dex.Perturbation(
    ...          [
    ...             dex.FluteMode(1e-5, lcfs, 2, 1, 0),
    ...             dex.FluteMode(1e-5, lcfs, 3, 1, 0),
    ...          ]
    ...     ),
    ... )
    >>>
    >>> strength = machine.quantity(50, "kiloVolt / meter")
    >>> mean = lcfs / 2
    >>> width = lcfs.value / 50
    >>>
    >>> # Now we can set the machine's electric field
    >>> machine.efield = dex.GaussianEfield(strength, mean, width)

    ```

    """

    _r: _PyEfield

    def __init__(
        self,
        strength: PlainQuantity[float] | float,
        mean: MagneticFlux,
        width: float,
    ) -> None:
        if isinstance(strength, PlainQuantity):
            strength = strength.to("NormVolt / NormMeter").m
        self._r = _PyEfield.build_gaussian_efield(
            strength=strength,
            mean=mean._r,
            width=width,
        )
        MachineObject.__init__(self)
        Efield.__init__(self)

    @property
    def strength(self) -> float:
        return self._r.strength

    @property
    def mean(self) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.mean)

    @property
    def width(self) -> float:
        return self._r.width


EfieldObject: TypeAlias = NoEfield | GaussianEfield
r"""Available [`Efield`][dexter.machine.base.Efield] objects.

+ [`NoEfield`][dexter.NoEfield]: No electric field.
+ [`GaussianEfield`][dexter.GaussianEfield]: Radial-Gaussian electric field.
"""
