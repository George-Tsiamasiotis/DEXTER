import pytest
import numpy as np
from math import isfinite
from math import pi as PI
from dexter import Perturbation, CosHarmonic, NcHarmonic


def test_invalid_perturbation(nc_harmonic: NcHarmonic):
    with pytest.raises(TypeError):
        Perturbation([CosHarmonic(1e-3, 1, 2, 0.0), nc_harmonic])  # type: ignore


def test_empty_perturbation():
    p = Perturbation([])
    assert p.harmonics == []
    assert len(p) == 0
    with pytest.raises(BaseException):
        p[0]

    args = (0.01, 1, 2, 3)
    assert p.p_of_psi(*args) == 0
    assert p.p_of_psip(*args) == 0
    assert p.dp_dpsi(*args) == 0
    assert p.dp_dpsip(*args) == 0
    assert p.dp_of_psi_dtheta(*args) == 0
    assert p.dp_of_psip_dzeta(*args) == 0
    assert p.dp_of_psi_dtheta(*args) == 0
    assert p.dp_of_psip_dzeta(*args) == 0
    assert p.dp_of_psi_dt(*args) == 0
    assert p.dp_of_psip_dt(*args) == 0


def test_cos_perturbation():
    p = Perturbation(
        [
            CosHarmonic(1e-3, 1, 2, 0.0),
            CosHarmonic(1e-3, 1, 3, 0.0),
            CosHarmonic(1e-3, 1, 4, 0.0),
            CosHarmonic(1e-3, 1, 5, 0.0),
        ]
    )
    assert isinstance(p.harmonics, list)
    assert len(p) == 4
    _test_perturbation_vectorized_evals(p)


def test_nc_perturbation(nc_perturbation: Perturbation):
    assert len(nc_perturbation) == 3
    assert isinstance(nc_perturbation.harmonics, list)
    _test_perturbation_vectorized_evals(nc_perturbation)


def _test_perturbation_vectorized_evals(perturbation: Perturbation):
    methods = [
        perturbation.p_of_psi,
        perturbation.p_of_psip,
        perturbation.dp_dpsi,
        perturbation.dp_dpsip,
        perturbation.dp_of_psi_dtheta,
        perturbation.dp_of_psip_dtheta,
        perturbation.dp_of_psi_dzeta,
        perturbation.dp_of_psip_dzeta,
        perturbation.dp_of_psi_dt,
        perturbation.dp_of_psip_dt,
    ]

    # 0D evaluations
    flux = 1e-5
    theta = PI
    zeta = 2 * PI
    t = 10
    for method in methods:
        assert isfinite(method(flux, theta, zeta, t))
        assert isinstance(method(flux, theta, zeta, t), float)

    # 1D Evaluations
    fluxes = np.linspace(1e-5, 1e-4, 5)
    thetas = np.linspace(0, PI, 5)
    zetas = np.linspace(0, 2 * PI, 5)
    ts = np.linspace(0, 100, 5)
    for method in methods:
        assert method(fluxes, thetas, zetas, ts).ndim == 1
        assert isinstance(method(fluxes, thetas, zetas, ts), np.ndarray)

    # 4D Evaluations
    flux_grid = np.random.random([2] * 4) * 1e-5
    theta_grid = np.random.random([2] * 4) * PI
    zeta_grid = np.random.random([2] * 4) * 2 * PI
    t_grid = np.random.random([2] * 4) * 1e-5
    assert flux_grid.ndim == 4
    for method in methods:
        assert method(flux_grid, theta_grid, zeta_grid, t_grid).ndim == 4
        assert isinstance(method(flux_grid, theta_grid, zeta_grid, t_grid), np.ndarray)
