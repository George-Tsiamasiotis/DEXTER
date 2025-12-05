from dexter import Qfactor


def test_fields(qfactor: Qfactor):
    assert isinstance(qfactor.path, str)
    assert isinstance(qfactor.typ, str)
    assert qfactor.psip_data.ndim == 1
    assert qfactor.q_data.ndim == 1
    assert qfactor.psi_data.ndim == 1


def test_eval(qfactor: Qfactor):
    psip = 0.015
    assert isinstance(qfactor.q(psip), float)
    assert isinstance(qfactor.psi(psip), float)
    assert isinstance(qfactor.dpsi_dpsip(psip), float)


def test_magic(qfactor: Qfactor):
    assert len(qfactor) == len(qfactor.psip_data)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)
