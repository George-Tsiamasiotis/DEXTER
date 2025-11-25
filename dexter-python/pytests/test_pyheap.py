import numpy as np
from dexter import Heap, HeapInitialConditions, Particle, MappingParameters
from dexter import Qfactor, Currents, Bfield, Perturbation


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
    assert heap.zetas.shape == (0, 0)
    assert heap.psips.shape == (0, 0)
    assert heap.thetas.shape == (0, 0)
    assert heap.psis.shape == (0, 0)
    str(heap)


def test_heap_poincare(
    qfactor: Qfactor,
    currents: Currents,
    bfield: Bfield,
    perturbation: Perturbation,
):
    psip = 0.09378

    num1 = 16
    psips = psip * np.ones(num1)
    zetas = np.linspace(-np.pi, np.pi, num1)

    initials = HeapInitialConditions(
        psips=psips,
        zetas=zetas,
        thetas=np.zeros(len(psips)),
        rhos=0.0001 * np.ones(len(psips)),
        mus=np.zeros(len(psips)),
    )

    heap = Heap(initials)

    params = MappingParameters("ConstTheta", 3.14, 10)

    heap.poincare(qfactor, currents, bfield, perturbation, params)
    print(heap)
