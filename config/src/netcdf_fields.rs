//! The names each variable is expected to appear in the netCDF file.
//!
//! If the naming convention changes, this is the only file we must update.

// ================== Scalars ==================

/// Magnetic field strength on the axis `B0` **in \[T\]**.
pub const BAXIS: &str = "baxis";
/// The horizontal position of the magnetic axis `R0` **in \[m\]**.
pub const RAXIS: &str = "raxis";
/// The vertical position of the magnetic axis **in \[m\]**.
pub const ZAXIS: &str = "zaxis";
/// The geometrical axis (device major radius) **in \[m\]**.
pub const RGEO: &str = "rgeo";

// ================= Coordinates =================

/// The boozer toroidal angle `θ` **in \[rads\]**.
pub const THETA: &str = "theta";
/// The poloidal flux `ψp` **in Normalized Units**.
pub const PSIP_NORM: &str = "psip_norm";
/// The toroidal flux `ψ` **in Normalized Units**.
pub const PSI_NORM: &str = "psi_norm";
/// The radial coordinate `r` **in Normalized Units**.
pub const R_NORM: &str = "r_norm";
/// The poloidal mode numbers `m`.
pub const M: &str = "m";
/// The toroidal mode numbers `n`.
pub const N: &str = "n";
/// The poloidal flux `ψp` **in \[Tm\]**.
pub const PSIP: &str = "psip";
/// The toroidal flux `ψ` **in \[Tm\]**.
pub const PSI: &str = "psi";
/// The radial coordinate r **in \[m\]**.
pub const R: &str = "r";

// ================ 1D Variables ================

/// q(ψp): The safety factor.
pub const Q: &str = "q";
/// g(ψp): The covariant toroidal plasma current **in \[Tm\]**.
pub const G: &str = "g";
/// I(ψp): The covariant poloidal plasma current **in \[Tm\]**.
pub const I: &str = "i";
/// g(ψp): The covariant toroidal plasma current **in Normalized Units**.
pub const G_NORM: &str = "g_norm";
/// I(ψp): The covariant poloidal plasma current **in Normalized Units**.
pub const I_NORM: &str = "i_norm";

// ================ 2D Variables ================

/// B(ψp, θ): The magnetic field strength in **in \[T\]**.
pub const B: &str = "b";
/// B(ψp, θ): The magnetic field strength in **in Normalized Units**.
pub const B_NORM: &str = "b_norm";
/// R(ψp, θ): The `R` coordinate with respect to boozer coordinates **in \[m\]**.
pub const RLAB: &str = "rlab";
/// Z(ψp, θ): The `Z` coordinate with respect to boozer coordinates **in \[m\]**.
pub const ZLAB: &str = "zlab";

// ================ 3D Variables ================

/// The 3D array containing all the `α{m,n}(ψp)` 1D arrays **in Normalized Units**.
pub const ALPHAS_NORM: &str = "alphas_norm";
/// The 3D array containing all the `α{m,n}(ψp)` 1D arrays **i \[m\]**.
pub const ALPHAS: &str = "alphas";
/// The 3D array containing all the `φ{m,n}(ψp)` 1D arrays.
pub const PHASES: &str = "phases";
