r""" """

from dexter.machine.machine import Machine

from dexter._utils import _ReprStrImpl, _RustTypeWrapper
from dexter._core import _PyStationaryCurveSegment, _PyStationaryCurve
from dexter.types import Array1, MagneticFluxKind


class StationaryCurveSegment(_ReprStrImpl, _RustTypeWrapper):
    r"""A single segment of the [`StationaryCurve`][dexter.StationaryCurve].

    Attributes
    ----------
    theta
        The $\theta$ array.
    flux
        The flux ($\psi$ or $\psi_p$) array.
    """

    _r = _PyStationaryCurveSegment

    @property
    def theta(self) -> Array1:
        return self._r.theta

    @property
    def flux(self) -> Array1:
        return self._r.flux

    def __len__(self) -> int:
        """Returns the total number of points."""
        return self._r.__len__()


class StationaryCurve(_ReprStrImpl):
    r"""The stationary curve of the Hamiltonian.

    The stationary curve is defined through the equation $\partial\mathcal{H}/\partial\theta=0$.
    In the absence of an electric field, this is equivalent to \partial B/\partial\theta=0$.

    The stationary curve is described by one or more segments of the form $f(\psi,\theta)=0$.

    Parameters
    ----------
    machine
        The machine under study.

    Attributes
    ----------
    flux_kind
        The kind of magnetic flux w.r.t. the stationary curve is defined.
    """

    _r: _PyStationaryCurve

    def __init__(self, machine: Machine) -> None:
        self._r = _PyStationaryCurve(
            machine.qfactor._r, machine.current._r, machine.bfield._r
        )

    def segments(self) -> list[StationaryCurveSegment]:
        r"""Returns a list of the comprising segments."""
        return [StationaryCurveSegment._wrap(segment) for segment in self._r.segments]

    @property
    def flux_kind(self) -> MagneticFluxKind:
        return self._r.flux_kind
