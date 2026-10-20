import numpy as np
import dexter as dex


def test_no_efield():
    efield = dex.NoEfield()
    _test_efield_base(efield)


def test_gaussian_efield(lar_machine: dex.Machine):
    strength = (
        lar_machine.quantity(750, "kiloVolt / meter").to("NormVolt / NormMeter").m
    )
    mean = lar_machine.psi_last / 2
    width = lar_machine.psi_last.value / 50
    efield = dex.GaussianEfield(strength, mean, width)

    assert efield.strength == strength
    assert efield.mean == mean
    assert efield.width == width

    _test_efield_base(efield)


def _test_efield_base(efield: dex.EfieldObject):

    efield.__repr__()
    efield.__str__()

    assert efield.machine_type in ["Numerical", "Analytical"]
    assert efield.psi_state in ["Good", "Bad"]
    assert efield.psip_state in ["Good", "Bad"]

    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)

    try:

        efield.eval_phi(psi=flux)
        efield.eval_phi(psi=fluxes)
        efield.eval_phi(psip=flux)
        efield.eval_phi(psip=fluxes)

        efield.eval_phi_deriv(psi=flux)
        efield.eval_phi_deriv(psi=fluxes)
        efield.eval_phi_deriv(psip=flux)
        efield.eval_phi_deriv(psip=fluxes)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )
