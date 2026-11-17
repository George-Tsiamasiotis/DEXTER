import pytest
import dexter as dex


def test_stepping_method():
    with pytest.raises(RuntimeError):
        dex.SteppingMethod()
    _energy = dex.SteppingMethod.EnergyAdaptiveStep(1e-4, 1e-10)
    _error = dex.SteppingMethod.ErrorAdaptiveStep(1e-4, 1e-10)
    _fixed = dex.SteppingMethod.FixedStep(0.01)


def test_solver_params():
    _default = dex.SolverParams()
    _energy = dex.SolverParams(method=dex.SteppingMethod.EnergyAdaptiveStep(1e-2, 1e-3))
    _error = dex.SolverParams(method=dex.SteppingMethod.ErrorAdaptiveStep(1e-2, 1e-3))
    _fixed = dex.SolverParams(method=dex.SteppingMethod.FixedStep(1e-2))

    dex.SolverParams(max_steps=10)
    dex.SolverParams(safety_factor=0.8)


def test_intersect_params():
    # fmt: off
    dex.IntersectParams(intersection="ConstZeta", angle=2, turns=1000, directionality="Initial")
    dex.IntersectParams(intersection="ConstZeta", angle=2, turns=1000, directionality="Both")
    dex.IntersectParams(intersection="ConstTheta", angle=2, turns=1000, directionality="DotNegative")
    dex.IntersectParams(intersection="ConstTheta", angle=2, turns=1000, directionality="DotPositive")
    # fmt: on
