import pytest
import dexter as dex


def test_stepping_method():
    with pytest.raises(RuntimeError):
        dex.SteppingMethod()
    energy = dex.SteppingMethod.EnergyAdaptiveStep(1e-4, 1e-10)
    error = dex.SteppingMethod.ErrorAdaptiveStep(1e-4, 1e-10)
    fixed = dex.SteppingMethod.FixedStep(0.01)


def test_solver_params():
    default = dex.SolverParams()
    energy = dex.SolverParams(dex.SteppingMethod.EnergyAdaptiveStep(1e-2, 1e-3))
    error = dex.SolverParams(dex.SteppingMethod.ErrorAdaptiveStep(1e-2, 1e-3))
    fixed = dex.SolverParams(dex.SteppingMethod.FixedStep(1e-2))

    dex.SolverParams(max_steps=10)
    dex.SolverParams(safety_factor=0.8)


def test_intersect_params():
    dex.IntersectParams("ConstZeta", 2, 1000, "Initial")
    dex.IntersectParams("ConstZeta", 2, 1000, "Both")
    dex.IntersectParams("ConstTheta", 2, 1000, "DotNegative")
    dex.IntersectParams("ConstTheta", 2, 1000, "DotPositive")
