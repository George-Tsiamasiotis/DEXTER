from dexter import NcQfactor, UnityQfactor


def test_nc_qfactor_getters(nc_qfactor: NcQfactor):
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.typ, str)
    assert nc_qfactor.psip_data.ndim == 1
    assert nc_qfactor.q_data.ndim == 1
    assert nc_qfactor.psi_data.ndim == 1
    assert len(nc_qfactor) == len(nc_qfactor.psip_data)
    assert isinstance(nc_qfactor.__str__(), str)
    assert isinstance(nc_qfactor.__repr__(), str)


def test_nc_qfactor_evals(nc_qfactor: NcQfactor):
    psip = 0.015
    assert isinstance(nc_qfactor.q(psip), float)
    assert isinstance(nc_qfactor.psi(psip), float)
    assert isinstance(nc_qfactor.dpsi_dpsip(psip), float)


def test_unity_qfactor(unity_qfactor: UnityQfactor):
    psip = 0.015
    assert unity_qfactor.q(psip) == 1
    assert unity_qfactor.psi(psip) == psip
    assert unity_qfactor.dpsi_dpsip(psip) == 1
    str(unity_qfactor)
