"""This file mirrors all the definitions made in the `py-dexter` Rust API.

Each class representing a dexter-machine trait also provides all possible getters of all
implementors. The higher level wrappers are responsible for re-exporting the correct ones. An
`InvalidVariant` Exception is raised when python accesses a field that the wrapped type does not
have. It should not be visible to the user unless the call `._r` explicitly.
"""

from numpy import nan as NAN

from dexter.types import (
    Array1,
    Array2,
    ArrayShape,
    EnergyPzetaPosition,
    MagneticFluxKind,
    FluxCoordinateState,
    IntegrationStatus,
    Interpolation1dType,
    Interpolation2dType,
    OrbitType,
    PhaseMethod,
    MachineType,
    CoordinateSet,
    SteppingMethod,
    Intersection,
)

class _PyMagneticFlux:

    value: float
    kind: MagneticFluxKind

    @classmethod
    def toroidal(cls, value: float) -> _PyMagneticFlux: ...
    @classmethod
    def poloidal(cls, value: float) -> _PyMagneticFlux: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyGeometry:  # Trait and all possible getters
    machine_type: MachineType
    psi_state: FluxCoordinateState
    psip_state: FluxCoordinateState
    baxis: float
    raxis: float
    zaxis: float
    rgeo: float
    rlast: float
    path: str
    netcdf_version: str
    interp1d_type: Interpolation1dType
    interp2d_type: Interpolation2dType
    rlab_last: Array1
    zlab_last: Array1
    shape: ArrayShape
    psi_last: _PyMagneticFlux | None
    psip_last: _PyMagneticFlux | None

    @classmethod
    def build_lar(cls, baxis: float, raxis: float, rlast: float) -> _PyGeometry: ...
    @classmethod
    def build_nc(
        cls,
        path: str,
        interp1d_type: Interpolation1dType,
        interp2d_type: Interpolation2dType,
    ) -> _PyGeometry: ...
    def eval_r(self, psi: float, psip: float) -> float: ...
    def eval_psi_of_r(self, r: float) -> float: ...
    def eval_psip_of_r(self, r: float) -> float: ...
    def eval_rlab(self, theta: float, psi: float, psip: float) -> float: ...
    def eval_zlab(self, theta: float, psi: float, psip: float) -> float: ...
    def eval_jacobian(self, theta: float, psi: float, psip: float) -> float: ...
    def get_array(self, name: str) -> Array1 | None: ...
    def get_array2d(self, name: str) -> Array2 | None: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyQfactor:  # Trait and all possible getters
    machine_type: MachineType
    psi_state: FluxCoordinateState
    psip_state: FluxCoordinateState
    psi_last: _PyMagneticFlux
    psip_last: _PyMagneticFlux
    qlast: float
    qaxis: float
    path: str
    netcdf_version: str
    interp_type: Interpolation1dType

    @classmethod
    def build_unity(cls, lcfs: _PyMagneticFlux) -> _PyQfactor: ...
    @classmethod
    def build_parabolic(
        cls,
        qaxis: float,
        qlast: float,
        lcfs: _PyMagneticFlux,
    ) -> _PyQfactor: ...
    @classmethod
    def build_nc(cls, path: str, interp_type: Interpolation1dType) -> _PyQfactor: ...
    def eval_q(self, psi: float, psip: float) -> float: ...
    def eval_other(self, psi: float, psip: float) -> float: ...
    def eval_psi_of_q(self, q: float) -> float: ...
    def eval_psip_of_q(self, q: float) -> float: ...
    def eval_deriv_of_other(self, psi: float, psip: float) -> float: ...
    def eval_deriv_wrt_other(self, psi: float, psip: float) -> float: ...
    def eval_iota(self, psi: float, psip: float) -> float: ...
    def get_array(self, name: str) -> Array1: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyCurrent:  # Trait and all possible getters
    machine_type: MachineType
    psi_state: FluxCoordinateState
    psip_state: FluxCoordinateState
    path: str
    netcdf_version: str
    interp_type: Interpolation1dType

    @classmethod
    def build_lar(cls) -> _PyCurrent: ...
    @classmethod
    def build_nc(cls, path: str, interp_type: Interpolation1dType) -> _PyCurrent: ...
    def eval_g(self, psi: float, psip: float) -> float: ...
    def eval_i(self, psi: float, psip: float) -> float: ...
    def eval_g_deriv(self, psi: float, psip: float) -> float: ...
    def eval_i_deriv(self, psi: float, psip: float) -> float: ...
    def get_array(self, name: str) -> Array1 | None: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyBfield:  # Trait and all possible getters
    machine_type: MachineType
    psi_state: FluxCoordinateState
    psip_state: FluxCoordinateState
    path: str
    netcdf_version: str
    interp_type: Interpolation2dType
    baxis: float
    padding: int
    padding_theta: float
    shape: ArrayShape
    shape_padded: ArrayShape

    @classmethod
    def build_lar(cls) -> _PyBfield: ...
    @classmethod
    def build_nc(
        cls,
        path: str,
        interp_type: Interpolation2dType,
        padding: int,
    ) -> _PyBfield: ...
    def eval_b(
        self,
        theta: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_flux(
        self,
        theta: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_theta(
        self,
        theta: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def get_array(self, name: str) -> Array1 | None: ...
    def get_array2d(self, name: str) -> Array2: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyMode:  # Trait and all possible getters
    machine_type: MachineType
    psi_state: FluxCoordinateState
    psip_state: FluxCoordinateState
    m: int
    n: int
    lcfs: _PyMagneticFlux
    epsilon: float
    phase: float
    path: str
    netcdf_version: str
    interp_type: Interpolation1dType
    phase_method: PhaseMethod
    analytical_threshold_index: int
    phase_average: float | None

    @classmethod
    def build_flute(
        cls,
        epsilon: float,
        lcfs: _PyMagneticFlux,
        m: int,
        n: int,
        phase: float,
    ) -> _PyMode: ...
    @classmethod
    def build_nc(
        cls,
        path: str,
        interp_type: Interpolation1dType,
        m: int,
        n: int,
        phase_method: PhaseMethod,
        analytical_threshold_index: int,
    ) -> _PyMode: ...
    def eval_amplitude(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_phase(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_m(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_flux(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_theta(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_zeta(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_t(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def get_array(self, name: str) -> Array1 | None: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PyPerturbation:
    def __init__(self, modes: list[_PyMode]) -> None: ...
    def eval_p(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_flux(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_theta(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_zeta(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def eval_deriv_t(
        self,
        theta: float,
        zeta: float,
        t: float,
        psi: float,
        psip: float,
    ) -> float: ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

# ================================================================================================

class _PyInitialConditions:
    t0: float
    flux0: _PyMagneticFlux
    theta0: float
    zeta0: float
    rho0: float | None
    pzeta0: float | None
    mu0: float
    coordinate_set: CoordinateSet

    @classmethod
    def boozer(
        cls,
        t0: float,
        flux0: _PyMagneticFlux,
        theta0: float,
        zeta0: float,
        rho0: float,
        mu0: float,
    ) -> _PyInitialConditions: ...
    @classmethod
    def mixed(
        cls,
        t0: float,
        flux0: _PyMagneticFlux,
        theta0: float,
        zeta0: float,
        pzeta0: float,
        mu0: float,
    ) -> _PyInitialConditions: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class _PySolverParams:
    def __init__(
        self,
        stepping_method: SteppingMethod | None,
        max_steps: int | None,
        first_step: float | None,
        safety_factor: float | None,
        energy_rel_tol: float | None,
        energy_abs_tol: float | None,
        error_rel_tol: float | None,
        error_abs_tol: float | None,
    ) -> None: ...

class _PyIntersectParams:
    def __init__(
        self,
        intersection: Intersection,
        angle: float,
        turns: int,
    ) -> None: ...

class _PyParticle:
    initial_conditions: _PyInitialConditions
    integration_status: IntegrationStatus
    steps_taken: int
    steps_stored: int
    duration: str
    initial_energy: float | None
    final_energy: float | None
    energy_var: float | None
    energy_pzeta_position: EnergyPzetaPosition
    orbit_type: OrbitType
    omega_theta: float | None
    omega_zeta: float | None
    qkinetic: float | None
    flux_cache_hits: int
    flux_cache_misses: int
    theta_cache_hits: int
    theta_cache_misses: int
    mode_cache_hits: int
    mode_cache_misses: int

    def __init__(self, initial: _PyInitialConditions) -> None: ...
    def integrate(
        self,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        teval: tuple[float, float],
        solver_params: _PySolverParams,
    ) -> None: ...
    def intersect(
        self,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        intersect_params: _PyIntersectParams,
        solver_params: _PySolverParams,
    ) -> None: ...
    def close(
        self,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
        perturbation: _PyPerturbation,
        periods: int,
        solver_params: _PySolverParams,
    ) -> None: ...
    def classify(
        self,
        qfactor: _PyQfactor,
        current: _PyCurrent,
        bfield: _PyBfield,
    ) -> None: ...
    def print_caches(self) -> None: ...
    def discard_arrays(self) -> None: ...
    def get_array(self, name: str) -> Array1: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
