import numpy as np
from dexter import Heap, HeapInitialConditions, Particle


def test_extraction():
    initials = HeapInitialConditions(
        thetas=np.zeros(5),
        psips=np.zeros(5),
        rhos=np.zeros(5),
        zetas=np.zeros(5),
        mus=np.zeros(5),
    )

    assert isinstance(initials.thetas, np.ndarray)
    assert isinstance(initials.psips, np.ndarray)
    assert isinstance(initials.rhos, np.ndarray)
    assert isinstance(initials.zetas, np.ndarray)
    assert isinstance(initials.mus, np.ndarray)
    assert len(initials) == 5


def test_heap_initialization():
    initials = HeapInitialConditions(
        thetas=np.zeros(5),
        psips=np.zeros(5),
        rhos=np.zeros(5),
        zetas=np.zeros(5),
        mus=np.zeros(5),
    )
    heap = Heap(initials)
    assert len(heap) == 5
    assert isinstance(heap[2], Particle)
