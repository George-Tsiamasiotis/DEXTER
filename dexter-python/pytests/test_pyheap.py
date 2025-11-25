import numpy as np
from dexter import HeapInitialConditions


def test_extraction():
    init = HeapInitialConditions(
        thetas=np.zeros(5),
        psips=np.zeros(5),
        rhos=np.zeros(5),
        zetas=np.zeros(5),
        mus=np.zeros(5),
    )

    assert isinstance(init.thetas, np.ndarray)
    assert isinstance(init.psips, np.ndarray)
    assert isinstance(init.rhos, np.ndarray)
    assert isinstance(init.zetas, np.ndarray)
    assert isinstance(init.mus, np.ndarray)
    assert len(init) == 5
