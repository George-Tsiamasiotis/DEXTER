import numpy as np

from dexter import Equilibrium, COMs, EnergyPzetaPlane, Parabola


def test_numerical_energy_pzeta_plane(nc_equilibrium: Equilibrium):

    coms = COMs(mu=1e-4)
    plane = coms.build_energy_pzeta_plane(nc_equilibrium)
    assert isinstance(plane, EnergyPzetaPlane)
    assert plane.tp_pzeta_interval.ndim == 1
    assert plane.tp_lower.ndim == 1
    assert plane.tp_upper.ndim == 1
    assert coms.mu == 1e-4
    _test_parabola(plane.axis_parabola)
    _test_parabola(plane.left_wall_parabola)
    _test_parabola(plane.right_wall_parabola)


def _test_parabola(parabola: Parabola):
    assert isinstance(parabola.a, float)
    assert isinstance(parabola.b, float)
    assert isinstance(parabola.c, float)
    assert isinstance(parabola.eval(1), float)

    x = np.linspace(-1, 1, 10)
    assert isinstance(parabola.eval_array1(x), np.ndarray)
