"""Defines the `StationaryCurve` and `StationaryCurveSegment` classes."""

from functools import cached_property

from dexter.types import Array1, MagneticFluxKind

from dexter.machine.machine import Machine
from dexter.machine.flux import MagneticFluxArray

from dexter._utils import _ReprStrImpl, _RustTypeWrapper
from dexter._core import _PyStationaryCurveSegment, _PyStationaryCurve


class StationaryCurve(_ReprStrImpl, _RustTypeWrapper):
    r"""The stationary curve of the Hamiltonian.

    The stationary curve is defined through the equation $\partial\mathcal{H}/\partial\theta = 0$.
    In the absence of an electric field, this is equivalent to \partial B/\partial\theta = 0$.

    The stationary curve is described by one or more parametrized segments of the form
    $\partial\mathcal{H}/\partial\theta = 0$

    Attributes
    ----------
    segments
        The segments comprising the stationary curve.
    """

    _r: _PyStationaryCurve

    def __init__(self, machine: Machine) -> None:
        self._r = _PyStationaryCurve(
            qfactor=machine.qfactor._r,
            current=machine.current._r,
            bfield=machine.bfield._r,
        )

    @cached_property  # immutable type
    def segments(self) -> list[StationaryCurveSegment]:
        """The segments comprising the stationary curve."""
        return [StationaryCurveSegment._wrap(_segment) for _segment in self._r.segments]

    def __len__(self) -> int:
        r"""Returns the number of [`segments`][dexter.StationaryCurveSegment] comprising the curve."""
        return self._r.__len__()


class StationaryCurveSegment(_ReprStrImpl, _RustTypeWrapper):
    r"""A 2D parametrized curve. Represents an individual segment of the Stationary Curve.

    The parameter $t$ is chosen to be the index of the corresponding $(\psi/\psi_p, \theta)$ point.

    Attributes
    ----------
    flux_kind
        The kind of the magnetic flux w.r.t. the segment is defined.
    flux_array
        The contained magnetic flux $\psi$ or $\psi_p$ array.
    theta_array
        The contained $\theta$ array.
    """

    _r: _PyStationaryCurveSegment

    def __init__(self) -> None:
        raise RuntimeError("'StationaryCurveSegment' cannot be constructed")

    @property
    def flux_kind(self) -> MagneticFluxKind:
        return self._r.flux_kind

    @cached_property  # immutable type
    def flux_array(self) -> MagneticFluxArray:
        return MagneticFluxArray._wrap(self._r.flux_array())

    @cached_property  # immutable type
    def theta_array(self) -> Array1:
        return self._r.theta_array()

    def __len__(self) -> int:
        r"""Returns the number of the contained $(\theta, \psi/\psi_p)$ points."""
        return self._r.__len__()
