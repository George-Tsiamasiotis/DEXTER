"""Defines the `StationaryCurve` and `ParametrizedSegment` classes."""

from typing import cast
from functools import cached_property

from dexter.types import Array1, CriticalPointKind, MagneticFluxKind

from dexter._utils import _ReprStrImpl, _RustTypeWrapper
from dexter._core import _PyParametrizedSegment, _PyStationaryCurve, _PyCriticalPoint

from dexter.machine.flux import MagneticFlux, MagneticFluxArray


class StationaryCurve(_ReprStrImpl, _RustTypeWrapper):
    r"""The stationary curve of the Hamiltonian.

    The stationary curve is defined through the equation $\partial\mathcal{H}/\partial\theta = 0$.
    In the absence of an electric field, this is equivalent to \partial B/\partial\theta = 0$.

    The stationary curve is described by one or more parametrized segments of the form
    $\partial\mathcal{H}/\partial\theta = 0$

    This type is built through
    [`GCHamiltonian.build_stationary_curve`][dexter.GCHamiltonian.build_stationary_curve].

    Attributes
    ----------
    segments
        The segments comprising the stationary curve.
    """

    _r: _PyStationaryCurve

    def __init__(self) -> None:
        raise RuntimeError(
            "'StationaryCurve' must be constructed through GCHamiltonian.build_stationary_curve"
        )

    @cached_property  # immutable type
    def segments(self) -> list[ParametrizedSegment]:
        """The segments comprising the stationary curve."""
        return [ParametrizedSegment._wrap(_segment) for _segment in self._r.segments]


class ParametrizedSegment(_ReprStrImpl, _RustTypeWrapper):
    r"""A 2D parametrized curve. Represents an individual segment of the Stationary Curve.

    The parameter $t$ is chosen to be the index of the corresponding $(\psi/\psi_p, \theta)$ point.

    Attributes
    ----------
    flux_kind
        The kind of the magnetic flux w.r.t. the segment is defined.
    t_array
        The parameter $t$ array.
    flux_array
        The contained magnetic flux $\psi$ or $\psi_p$ array.
    theta_array
        The contained $\theta$ array.
    theta_dot_array
        The $\dot\theta$ array, if it has been calculated
    """

    _r: _PyParametrizedSegment

    def __init__(self) -> None:
        raise RuntimeError("'ParametrizedSegment' cannot be constructed")

    @property
    def flux_kind(self) -> MagneticFluxKind:
        return self._r.flux_kind

    @cached_property  # immutable type
    def flux_array(self) -> MagneticFluxArray:
        # always exists
        return MagneticFluxArray._wrap(self._r.flux_array())

    @cached_property  # immutable type
    def t_array(self) -> Array1:
        # always exists
        return cast(Array1, self._r.get_array("t_array"))

    @cached_property  # immutable type
    def theta_array(self) -> Array1:
        # always exists
        return cast(Array1, self._r.get_array("theta_array"))

    @property  # immutable type
    def theta_dot_array(self) -> Array1 | None:
        return self._r.get_array("theta_dot_array")

    def __len__(self) -> int:
        r"""Returns the number of the contained $(\theta, \psi/\psi_p)$ points."""
        return self._r.__len__()


class CriticalPoint(_ReprStrImpl, _RustTypeWrapper):
    r"""A critical point of the Hamiltonian.

    A simple container type that is created through
    [`EnergyPzetaPlane.calculate_critical_points`][dexter.EnergyPzetaPlane.calculate_critical_points].

    Attributes
    ----------
    theta
        The point's $\theta$ position.
    flux
        The point's $\psi/\psi_p$ position.
    energy
        The point's energy.
    pzeta
        The point's $P_\zeta$.
    mu
        The point's $\mu$.
    kind
        The point's kind.
    """

    _r: _PyCriticalPoint

    @property
    def theta(self) -> float:
        return self._r.theta

    @property
    def flux(self) -> MagneticFlux:
        return MagneticFlux._wrap(self._r.flux)

    @property
    def pzeta(self) -> float:
        return self._r.pzeta

    @property
    def energy(self) -> float:
        return self._r.energy

    @property
    def mu(self) -> float:
        return self._r.mu

    @property
    def kind(self) -> CriticalPointKind:
        return self._r.kind
