"""This file mirrors all the definitions made in the `py-dexter` Rust API."""

from typing import Optional, TypeAlias

from dexter.types import (
    ArrayShape,
    Array1,
    Array2,
    NetCDFVersion,
    EquilibriumType,
    Interp1DType,
    Interp2DType,
    FluxCoordinate,
    FluxState,
    PhaseMethod,
    CoordinateSet,
    Intersection,
    IntegrationStatus,
    OrbitType,
    SteppingMethod,
    Routine,
)

# ================================================================================================

def _py_get_max_threads() -> int:
    """PyO3 export of `dexter_common::get_max_threads`."""

def _py_set_num_threads(num: int):
    """PyO3 export of `dexter_common::set_num_threads`."""

# ================================================================================================

class _PyLastClosedFluxSurface:
    """PyO3 export of `py_dexter::PyLastClosedFluxSurface`."""

    kind: FluxCoordinate
    value: float

    def __init__(self, kind: FluxCoordinate, value: float): ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyLarGeometry:
    """PyO3 export of `dexter_equilibrium::LarGeometry`."""

    equilibrium_type: EquilibriumType
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rlast: float
    psi_last: float
    psi_state: FluxState
    psip_state: FluxState
    rlab_last: Array1
    zlab_last: Array1

    def __init__(self, baxis: float, raxis: float, rlast: float) -> None: ...
    def r_of_psi(self, psi: float) -> float: ...
    def r_of_psip(self, psip: float) -> float: ...
    def psi_of_r(self, r: float) -> float: ...
    def psip_of_r(self, r: float) -> float: ...
    def rlab_of_psi(self, psi: float, theta: float) -> float: ...
    def rlab_of_psip(self, psip: float, theta: float) -> float: ...
    def zlab_of_psi(self, psi: float, theta: float) -> float: ...
    def zlab_of_psip(self, psip: float, theta: float) -> float: ...
    def jacobian_of_psi(self, psi: float, theta: float) -> float: ...
    def jacobian_of_psip(self, psip: float, theta: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcGeometry:
    """PyO3 export of `dexter_equilibrium::NcGeometry`."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp1d_type: Interp1DType
    interp2d_type: Interp2DType
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rlast: float
    shape: ArrayShape
    psi_last: float
    psip_last: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    r_array: Array1
    rlab_array: Array2
    zlab_array: Array2
    jacobian_array: Array2
    rlab_last: Array1
    zlab_last: Array1

    def __init__(
        self,
        path: str,
        interp1d_type: Interp1DType,
        interp2d_type: Interp2DType,
    ) -> None: ...
    def psip_of_psi(self, psi: float) -> float: ...
    def psi_of_psip(self, psip: float) -> float: ...
    def r_of_psi(self, psi: float) -> float: ...
    def r_of_psip(self, psip: float) -> float: ...
    def psi_of_r(self, r: float) -> float: ...
    def psip_of_r(self, r: float) -> float: ...
    def rlab_of_psi(self, psi: float, theta: float) -> float: ...
    def rlab_of_psip(self, psip: float, theta: float) -> float: ...
    def zlab_of_psi(self, psi: float, theta: float) -> float: ...
    def zlab_of_psip(self, psip: float, theta: float) -> float: ...
    def jacobian_of_psi(self, psi: float, theta: float) -> float: ...
    def jacobian_of_psip(self, psip: float, theta: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyUnityQfactor:
    """PyO3 export of `dexter_equilibrium::UnityQfactor`."""

    equilibrium_type: EquilibriumType
    psi_last: float
    psip_last: float
    qlast: float
    qaxis: float
    psi_state: FluxState
    psip_state: FluxState

    def __init__(self, lcfs: _PyLastClosedFluxSurface): ...
    def psip_of_psi(self, psi: float) -> float: ...
    def psi_of_psip(self, psip: float) -> float: ...
    def q_of_psi(self, psi: float) -> float: ...
    def q_of_psip(self, psip: float) -> float: ...
    def dpsi_dpsip(self, psip: float) -> float: ...
    def dpsip_dpsi(self, psi: float) -> float: ...
    def iota_of_psi(self, psi: float) -> float: ...
    def iota_of_psip(self, psip: float) -> float: ...
    def psi_of_q(self, q: float) -> float: ...
    def psip_of_q(self, q: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyParabolicQfactor:
    """PyO3 export of `dexter_equilibrium::ParabolicQfactor`."""

    equilibrium_type: EquilibriumType
    psi_last: float
    psip_last: float
    qlast: float
    qaxis: float
    psi_state: FluxState
    psip_state: FluxState

    def __init__(
        self,
        qaxis: float,
        qlast: float,
        lcfs: _PyLastClosedFluxSurface,
    ) -> None: ...
    def psip_of_psi(self, psi: float) -> float: ...
    def psi_of_psip(self, psip: float) -> float: ...
    def q_of_psi(self, psi: float) -> float: ...
    def q_of_psip(self, psip: float) -> float: ...
    def dpsi_dpsip(self, psip: float) -> float: ...
    def dpsip_dpsi(self, psi: float) -> float: ...
    def iota_of_psi(self, psi: float) -> float: ...
    def iota_of_psip(self, psip: float) -> float: ...
    def psi_of_q(self, q: float) -> float: ...
    def psip_of_q(self, q: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcQfactor:
    """PyO3 export of `dexter_equilibrium::NcQfactor`."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    qaxis: float
    qlast: float
    psi_last: float
    psip_last: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    q_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None: ...
    def psip_of_psi(self, psi: float) -> float: ...
    def psi_of_psip(self, psip: float) -> float: ...
    def q_of_psi(self, psi: float) -> float: ...
    def q_of_psip(self, psip: float) -> float: ...
    def dpsi_dpsip(self, psip: float) -> float: ...
    def dpsip_dpsi(self, psi: float) -> float: ...
    def iota_of_psi(self, psi: float) -> float: ...
    def iota_of_psip(self, psip: float) -> float: ...
    def psi_of_q(self, q: float) -> float: ...
    def psip_of_q(self, q: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyLarCurrent:
    """PyO3 export of `dexter_equilibrium::LarCurrent`."""

    equilibrium_type: EquilibriumType
    psi_state: FluxState
    psip_state: FluxState

    def __init__(self) -> None: ...
    def g_of_psi(self, psi: float) -> float: ...
    def g_of_psip(self, psip: float) -> float: ...
    def i_of_psi(self, psi: float) -> float: ...
    def i_of_psip(self, psip: float) -> float: ...
    def dg_dpsi(self, psi: float) -> float: ...
    def dg_dpsip(self, psip: float) -> float: ...
    def di_dpsi(self, psi: float) -> float: ...
    def di_dpsip(self, psip: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcCurrent:
    """PyO3 export of `dexter_equilibrium::NcCurrent`."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    psi_last: float
    psip_last: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def __init__(self, path: str, interp_type: Interp1DType) -> None: ...
    def g_of_psi(self, psi: float) -> float: ...
    def g_of_psip(self, psip: float) -> float: ...
    def i_of_psi(self, psi: float) -> float: ...
    def i_of_psip(self, psip: float) -> float: ...
    def dg_dpsi(self, psi: float) -> float: ...
    def dg_dpsip(self, psip: float) -> float: ...
    def di_dpsi(self, psi: float) -> float: ...
    def di_dpsip(self, psip: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyLarBfield:
    """PyO3 export of `dexter_equilibrium::LarBfield`."""

    equilibrium_type: EquilibriumType
    psi_state: FluxState
    psip_state: FluxState

    def __init__(self) -> None: ...
    def b_of_psi(self, psi: float, theta: float) -> float: ...
    def b_of_psip(self, psip: float, theta: float) -> float: ...
    def db_dpsi(self, psi: float, theta: float) -> float: ...
    def db_dpsip(self, psip: float, theta: float) -> float: ...
    def db_of_psi_dtheta(self, psi: float, theta: float) -> float: ...
    def db_of_psip_dtheta(self, psip: float, theta: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcBfield:
    """PyO3 export of `dexter_equilibrium::NcBfield`."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    baxis: float
    padding: int
    shape: ArrayShape
    shape_padded: ArrayShape
    psi_last: float
    psip_last: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    b_array: Array2
    theta_array_padded: Array1
    b_array_padded: Array2

    def __init__(
        self,
        path: str,
        interp_type: Interp2DType,
        *,
        padding: int = 15,
    ) -> None: ...
    def b_of_psi(self, psi: float, theta: float) -> float: ...
    def b_of_psip(self, psip: float, theta: float) -> float: ...
    def db_dpsi(self, psi: float, theta: float) -> float: ...
    def db_dpsip(self, psip: float, theta: float) -> float: ...
    def db_of_psi_dtheta(self, psi: float, theta: float) -> float: ...
    def db_of_psip_dtheta(self, psip: float, theta: float) -> float: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyCosHarmonic:
    """PyO3 export of `dexter_equilibrium::CosHarmonic`."""

    equilibrium_type: EquilibriumType
    epsilon: float
    lcfs: _PyLastClosedFluxSurface
    phase: float
    m: int
    n: int
    psi_last: Optional[float]
    psip_last: Optional[float]
    psi_state: FluxState
    psip_state: FluxState

    def __init__(
        self,
        epsilon: float,
        lcfs: _PyLastClosedFluxSurface,
        m: int,
        n: int,
        phase: float,
    ) -> None: ...
    # fmt: off
    def alpha_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def alpha_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def phase_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def phase_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def h_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def h_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dtheta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dtheta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dzeta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    # fmt: on
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcHarmonic:
    """PyO3 export of `dexter_equilibrium::NcHarmonic`."""

    path: str
    netcdf_version: NetCDFVersion
    equilibrium_type: EquilibriumType
    interp_type: Interp1DType
    m: int
    n: int
    phase_method: PhaseMethod
    phase_average: float
    psi_phase_resonance: float
    psip_phase_resonance: float
    psi_last: float
    psip_last: float
    psi_state: FluxState
    psip_state: FluxState
    psi_array: Array1
    psip_array: Array1
    alpha_array: Array1
    phase_array: Array1

    def __init__(
        self,
        path: str,
        interp_type: Interp1DType,
        m: int,
        n: int,
        phase_method: PhaseMethod = "Zero",
    ) -> None: ...
    # fmt: off
    def alpha_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def alpha_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def phase_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def phase_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def h_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def h_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dtheta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dtheta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dzeta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dh_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    # fmt: on
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyCosPerturbation:
    """PyO3 export of `dexter_equilibrium::Perturbation<CosHarmonic>`."""

    harmonics: list[_PyCosHarmonic]

    def __init__(self, harmonics: list[_PyCosHarmonic]) -> None: ...
    # fmt: off
    def p_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def p_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dtheta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dtheta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dzeta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def __getitem__(self, index: int): ...
    def __len__(self) -> int: ...
    # fmt: on
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyNcPerturbation:
    """PyO3 export of `dexter_equilibrium::Perturbation<NcHarmonic>`."""

    harmonics: list[_PyNcHarmonic]

    def __init__(self, harmonics: list[_PyNcHarmonic]) -> None: ...
    # fmt: off
    def p_of_psi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def p_of_psip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_dpsi(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_dpsip(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dtheta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dtheta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dzeta(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dzeta(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psi_dt(self, psi: float, theta: float, zeta: float, t: float) -> float: ...
    def dp_of_psip_dt(self, psip: float, theta: float, zeta: float, t: float) -> float: ...
    # fmt: on
    def __getitem__(self, index: int): ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================
# ================================================================================================

class _PyParabola:
    """PyO3 export of `parabola::Parabola`.

    Only methods that are needed are exported.
    """

    a: float
    b: float
    c: float

    def __init__(self) -> None: ...
    def eval(self, x: float) -> float: ...
    def eval_array1(self, arr: Array1) -> Array1: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyEnergyPzetaPlane:
    """PyO3 export of `dexter_simulate::EnergyPzetaPlane`."""

    axis_parabola: _PyParabola
    left_wall_parabola: _PyParabola
    right_wall_parabola: _PyParabola
    tp_pzeta_interval: Array1
    tp_upper: Array1
    tp_lower: Array1
    mu: float

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyCOMs:
    """PyO3 export of `dexter_simulate::COMs`."""

    energy: Optional[float]
    pzeta: Optional[float]
    mu: Optional[float]

    def __init__(
        self,
        energy: Optional[float] = None,
        pzeta: Optional[float] = None,
        mu: Optional[float] = None,
    ) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyInitialFlux:
    """PyO3 export of `dexter_simulate::InitialFlux`."""

    kind: FluxCoordinate
    value: float

    def __init__(self, kind: FluxCoordinate, value: float) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyInitialConditions:
    """PyO3 export of `dexter_simulate::InitialConditions`."""

    t0: float
    flux0: _PyInitialFlux
    theta0: float
    zeta0: float
    mu0: float

    rho0: Optional[float]
    pzeta0: Optional[float]

    coordinate_set: CoordinateSet

    def __init__(self) -> None:
        raise RuntimeError("This object should not be constructed directly")

    @classmethod
    def boozer(
        cls,
        t0: float,
        flux0: _PyInitialFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ) -> _PyInitialConditions: ...
    @classmethod
    def mixed(
        cls,
        t0: float,
        flux0: _PyInitialFlux,
        theta0: float,
        zeta0: float,
        pzeta0: float,
        mu0: float,
    ) -> _PyInitialConditions: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyIntersectParams:
    """PyO3 export of `dexter_simulate::IntersectParams`."""

    intersection: Intersection
    angle: float
    turns: int

    def __init__(
        self,
        intersection: Intersection,
        angle: float,
        turns: int,
    ) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

_PyQfactor: TypeAlias = _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor
_PyCurrent: TypeAlias = _PyLarCurrent | _PyNcCurrent
_PyBfield: TypeAlias = _PyLarBfield | _PyNcBfield
_PyPerturbation: TypeAlias = _PyCosPerturbation | _PyNcPerturbation

class _PyParticle:
    """PyO3 export of `dexter_simulate::Particle`."""

    initial_conditions: _PyInitialConditions
    integration_status: IntegrationStatus
    steps_taken: int
    steps_stored: int
    initial_energy: Optional[float]
    final_energy: Optional[float]
    energy_var: Optional[float]
    orbit_type: OrbitType
    omega_theta: Optional[float]
    omega_zeta: Optional[float]
    qkinetic: Optional[float]
    t_array: Array1
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    zeta_array: Array1
    rho_array: Array1
    mu_array: Array1
    ptheta_array: Array1
    pzeta_array: Array1
    energy_array: Array1

    def __init__(self, initial_conditions: _PyInitialConditions) -> None: ...
    def integrate(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        teval: tuple[float, float],
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def intersect(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        intersect_params: _PyIntersectParams,
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def close(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        periods: int,
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

# ================================================================================================

class _PyInitialFluxArray:
    """PyO3 export of `py-dexter::PyInitialFluxArray`.

    Helper object to initialize an array of InitialFluxes.
    """

    kind: FluxCoordinate
    values: Array1

    def __init__(self, kind: FluxCoordinate, values: Array1): ...

class _PyQueueInitialConditions:
    """PyO3 export of `dexter_simulate::QueueInitialConditions`."""

    t_array: Array1
    flux_array: Array1
    theta_array: Array1
    zeta_array: Array1
    mu_array: Array1

    rho_array: Optional[Array1]
    pzeta_array: Optional[Array1]

    def __init__(self) -> None:
        raise RuntimeError("This object should not be constructed directly")

    @classmethod
    def boozer(
        cls,
        t0: Array1,
        flux0: _PyInitialFluxArray,
        theta0: Array1,
        zeta0: Array1,
        rho0: Array1,
        mu0: Array1,
    ) -> _PyQueueInitialConditions: ...
    @classmethod
    def mixed(
        cls,
        t0: Array1,
        flux0: _PyInitialFluxArray,
        theta0: Array1,
        zeta0: Array1,
        pzeta0: Array1,
        mu0: Array1,
    ) -> _PyQueueInitialConditions: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> _PyInitialConditions: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class _PyQueue:
    """PyO3 export of `dexter_simulate::Queue`."""

    initial_conditions: _PyQueueInitialConditions
    particle_count: int
    particles: list[_PyParticle]
    routine: Routine
    energy_array: Array1
    omega_theta_array: Array1
    omega_zeta_array: Array1
    qkinetic_array: Array1
    steps_taken_array: Array1
    steps_stored_array: Array1
    _durations_as_nanos: Array1

    def __init__(self, initial_conditions: _PyQueueInitialConditions) -> None: ...
    @classmethod
    def from_particles(cls, particles: list[_PyParticle]) -> _PyQueue: ...
    def integrate(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        teval: tuple[float, float],
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def intersect(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        intersect_params: _PyIntersectParams,
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def close(
        self,
        /,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        periods: int,
        *,
        stepping_method: Optional[SteppingMethod],
        max_steps: Optional[int],
        first_step: Optional[float],
        safety_factor: Optional[float],
        energy_rel_tol: Optional[float],
        energy_abs_tol: Optional[float],
        error_rel_tol: Optional[float],
        error_abs_tol: Optional[float],
    ): ...
    def __getitem__(self, index: int) -> _PyParticle: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
