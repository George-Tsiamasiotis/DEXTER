from dexter import Harmonic


def test_derived_fields(harmonic1: Harmonic):
    assert isinstance(harmonic1.path, str)
    assert isinstance(harmonic1.typ, str)
    assert isinstance(harmonic1.phase_average, float)
    assert isinstance(harmonic1.m, int)
    assert isinstance(harmonic1.n, int)
    assert harmonic1.psip_data.ndim == 1
    assert harmonic1.a_data.ndim == 1
    assert harmonic1.phase_data.ndim == 1


def test_eval(harmonic1: Harmonic):
    (psip, theta, zeta) = 0.015, 3.14, 0
    assert isinstance(harmonic1.h(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dpsip(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dtheta(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dzeta(psip, theta, zeta), float)
    assert isinstance(harmonic1.dh_dt(psip, theta, zeta), float)
    assert isinstance(harmonic1.a(psip), float)
    assert isinstance(harmonic1.da_dpsip(psip), float)
    assert isinstance(harmonic1.phase(psip), float)


def test_magic(harmonic1: Harmonic):
    assert len(harmonic1) == len(harmonic1.psip_data)
    assert isinstance(harmonic1.__str__(), str)
    assert isinstance(harmonic1.__repr__(), str)
