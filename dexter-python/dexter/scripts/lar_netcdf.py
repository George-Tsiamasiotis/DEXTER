"""Creates a stub .npz file, with data produces from the analytical LAR formulas.

A parabolic curve is used for the q-factor profile.

Also adds 2 trivial perturbation

This is only intended for testing the `equilibrium` crate.
"""

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "output",
    help="the output 'npz' file",
)
args = parser.parse_args()

# ==========================================================================

import numpy as np
import xarray as xr
from pathlib import Path

OUTPUT = Path(args.output)

FLUX_SURFACES = 100
THETAS = 200

#########
# Scalars
#########

baxis = 1.5
raxis = 1.75
zaxis = 0
rgeo = 1.7

###################
# Parabolic Qfactor
###################

q0 = 1.1
qwall = 1.9
a = 0.5  # minor radius in meters
psi_wall = (a / raxis) ** 2 / 2


def parabolic_q(psi: np.ndarray) -> np.ndarray:
    return q0 + (qwall - q0) * (psi / psi_wall) ** 2


def parabolic_psip(psi: np.ndarray) -> np.ndarray:
    num = psi_wall / np.sqrt(q0 * (qwall - q0))
    atan = np.atan(psi * np.sqrt(qwall - q0) / (psi_wall * np.sqrt(q0)))
    return num * atan


psi_norm = np.linspace(0, psi_wall, FLUX_SURFACES)
q = parabolic_q(psi_norm)
psip_norm = parabolic_psip(psi_norm)
r_norm = np.sqrt(2 * psi_norm)

r = np.sqrt(2 * psi_norm) * raxis
psip = psip_norm * (baxis * raxis**2)
psi = psi_norm * (baxis * raxis**2)

############
# LAR Bfield
############

theta = np.linspace(0, 2 * np.pi, THETAS)


def lar_bfield(psi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    return 1 - np.sqrt(2 * psi) * np.cos(theta)


psi_norm_grid, theta_grid = np.meshgrid(psi_norm, theta)
b_norm = lar_bfield(psi_norm_grid, theta_grid)
b = b_norm / baxis

##########
# Currents
##########

g_norm = np.ones(FLUX_SURFACES)
i_norm = np.zeros(FLUX_SURFACES)
g = g_norm * baxis * raxis
i = i_norm * baxis * raxis

###############
# Perturbations
###############

m = np.asarray([1, 2])
n = np.asarray([2, 3])
# Should be just small variations
a12 = a13 = a22 = a23 = 1e-3 * np.cos(psip_norm)
p12 = p13 = p22 = p23 = np.sin(psip_norm)
alphas_norm = np.asarray([[a12, a23], [a13, a23]])
phases = np.asarray([[p12, p23], [p13, p23]])

alphas = alphas_norm * raxis

# Inject the (1, 2) mode with specific values to test that the extraction
# is done correctly. If either the indexes or the values change, the test
# must be updated accordingly.
alphas[0, 1, 0] = 1111
alphas[0, 1, -1] = 11111
phases[0, 1, 0] = 9999
phases[0, 1, -1] = 99999
alphas_norm[0, 1, 0] = 2222
alphas_norm[0, 1, -1] = 22222

#################
# Lab coordinates
#################

psip_norm_grid, theta_grid = np.meshgrid(psip_norm, theta)

rlab = rgeo + psi_norm_grid * np.cos(theta_grid)
zlab = psi_norm_grid * np.sin(theta_grid)


# ======================= Scalars


baxis_var = xr.Variable(
    data=baxis,
    dims=[],
    attrs=dict(
        description="The magnetic field strength on the axis `B0`",
        units="[T]",
        normalization="This value should be used for normalizations",
    ),
)

raxis_var = xr.Variable(
    data=raxis,
    dims=[],
    attrs=dict(
        description="The horizontal position of the magnetic axis R0",
        units="[m]",
        normalization="This value should be used for normalizations",
    ),
)

zaxis_var = xr.Variable(
    data=zaxis,
    dims=[],
    attrs=dict(description="The vertical position of the magnetic axis", units="[m]"),
)

rgeo_var = xr.Variable(
    data=rgeo,
    dims=[],
    attrs=dict(description="The geometrical axis (device major radius)", units="[m]"),
)


# ======================= Coordinates


theta_coord = xr.Variable(
    data=theta,
    dims=["theta"],
    attrs=dict(description="Boozer theta coordinate", units="[rads]"),
)

psip_norm_coord = xr.Variable(
    data=psip_norm,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal flux coordinate", units="Normalized"),
)

psi_norm_coord = xr.Variable(
    data=psi_norm,
    dims=["psi_norm"],
    attrs=dict(description="Toroidal flux coordinate", units="Normalized"),
)

m_coord = xr.Variable(
    data=m,
    dims=["m"],
    attrs=dict(description="Poloidal mode number"),
)

n_coord = xr.Variable(
    data=n,
    dims=["n"],
    attrs=dict(description="Toroidal mode number"),
)

r_norm_coord = xr.Variable(
    data=r_norm,
    dims=["r_norm"],
    attrs=dict(description="Radial coordinate", units="Normalized"),
)


# ======================= SI Variables

psip_var = xr.Variable(
    data=psip,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal flux coordinate", units="[Tm]"),
)

r_var = xr.Variable(
    data=r,
    dims=["r_norm"],
    attrs=dict(description="Radial coordinate", units="[m]"),
)

q_var = xr.Variable(
    data=q,
    dims=["psip_norm"],
    attrs=dict(description="Magnetic q-factor", units="dimensionless"),
)

psi_var = xr.Variable(
    data=psi,
    dims=["psi_norm"],
    attrs=dict(description="Toroidal flux coordinate", units="[Tm]"),
)


g_var = xr.Variable(
    data=g,
    dims=["psip_norm"],
    attrs=dict(description="Toroidal current", units="[Tm]"),
)

i_var = xr.Variable(
    data=i,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal current", units="[Tm]"),
)

b_var = xr.Variable(
    data=b.T,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="[T]"),
)

jacobian_var = xr.Variable(
    data=b_norm.T,  # FIXME:
    dims=["psip_norm", "theta"],
    attrs=dict(description="VMEC to Boozer Jacobian", units="[m/T]"),
)

rlab_var = xr.Variable(
    data=rlab.T,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Lab horizontal coordinate", units="[m]"),
)

zlab_var = xr.Variable(
    data=zlab.T,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Lab vertical coordinate", units="[m]"),
)


# ======================= Normalized Variables


g_norm_var = xr.Variable(
    data=g_norm,
    dims=["psip_norm"],
    attrs=dict(description="Toroidal current normalized", units="Normalized"),
)

i_norm_var = xr.Variable(
    data=i_norm,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal current normalized", units="Normalized"),
)

b_norm_var = xr.Variable(
    data=b_norm.T,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="Normalized"),
)


# ======================= Perturbations


alphas_var = xr.Variable(
    data=alphas,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="[m]"),
)


phases_var = xr.Variable(
    data=phases,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode phases φ{m,n}(ψp)", units="[rads]"),
)

alphas_norm_var = xr.Variable(
    data=alphas_norm,
    dims=["m", "n", "psip_norm"],
    attrs=dict(description="Mode amplitudes α{m,n}(ψp)", units="Normalized"),
)


# ========================================================


COORDS = {
    "psip_norm": psip_norm_coord,
    "theta": theta_coord,
    "m": m_coord,
    "n": n_coord,
    "psi_norm": psi_norm_coord,
    "r_norm": r_norm,
}

VARIABLES = {
    "baxis": baxis_var,
    "raxis": raxis_var,
    "zaxis": zaxis_var,
    "rgeo": rgeo_var,
    "psip": psip_var,
    "psi": psi_var,
    "r": r_var,
    "q": q_var,
    "g": g_var,
    "i": i_var,
    "b": b_var,
    "jacobian": jacobian_var,
    "rlab": rlab_var,
    "zlab": zlab_var,
    "alphas": alphas_var,
    "phases": phases_var,
    "g_norm": g_norm_var,
    "i_norm": i_norm_var,
    "b_norm": b_norm_var,
    "alphas_norm": alphas_norm_var,
}

dataset = xr.Dataset(
    data_vars=VARIABLES,
    coords=COORDS,
)

print("Created dataset")

dataset.close()
dataset.to_netcdf(OUTPUT)

print(f"Stored dataset at '{OUTPUT.absolute()}'")

raise SystemExit
