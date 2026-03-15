# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "semver",
#     "numpy",
#     "xarray",
# ]
# ///

"""Converts a Falcon netCDF to a DEXTER netCDF."""

import argparse
from typing import Any
import numpy as np
import xarray as xr
import tomllib
from semver import Version
from pathlib import Path
from datetime import datetime

with open("./Cargo.toml", "rb") as f:
    VERSION = Version.parse(tomllib.load(f)["workspace"]["metadata"]["netcdf_version"])

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=f"netCDF version {VERSION}",
)
parser.add_argument(
    "input",
    help="the input file name",
    type=str,
)
parser.add_argument(
    "output",
    help="the output file name",
    type=str,
)
args = parser.parse_args()

# ==========================================================================

INPUT = Path(args.input)
OUTPUT = Path(args.output)

BAXIS = "Baxis"
RAXIS = "raxis"
ZAXIS = "zaxis"
RGEO = ""
THETA = "boozer_theta"
PSI_NORM = "psi"
M = ""
N = ""
PSIP_NORM = ""
R_NORM = ""

PSIP = ""
R = ""
Q = "q"
PSI = "psi"
G = "g"
I = "i"
B = "b"
JACOBIAN = ""
RLAB = "R"
ZLAB = "Z"
G_NORM = "g_norm"
I_NORM = "I_norm"
B_NORM = "b_field_norm"

ALPHAS = ""
PHASES = ""
ALPHAS_NORM = ""

falcon = xr.open_dataset(INPUT)
shape = getattr(falcon, B_NORM).data.shape


def get(field: str, default: Any):
    raw = getattr(falcon, field, None)
    if raw is not None:
        return raw.data
    else:
        return default


def nan_array(shape: tuple[int, ...]) -> np.ndarray:
    return np.full(shape, np.nan)


baxis = get(BAXIS, None)
raxis = get(RAXIS, None)
zaxis = get(ZAXIS, None)
rgeo = get(RGEO, None)

theta = get(THETA, nan_array(shape[1]))
psi_norm = np.insert(get(PSI_NORM, nan_array(shape[0])), 0, 0)
m = np.asarray([])
n = np.asarray([])
psip_norm = np.insert(get(PSIP_NORM, nan_array(shape[0])), 0, 0)
r_norm = np.insert(get(R_NORM, nan_array(shape[0])), 0, 0)

psip = np.insert(get(PSIP, nan_array(shape[0])), 0, 0)
r = np.insert(get(R, nan_array(shape[0])), 0, 0)
_q = get(Q, nan_array(shape[0]))
q = np.insert(_q, 0, _q[0])
psi = np.insert(get(PSI, nan_array(shape[0])), 0, 0)
_g = get(G, nan_array(shape[0]))
g = np.insert(_g, 0, _g[0])
i = np.insert(get(I, nan_array(shape[0])), 0, 0)
b = np.insert(get(B, nan_array(shape)), 0, np.full(shape[1], baxis), axis=0)
rlab = np.insert(get(RLAB, nan_array(shape)), 0, np.full(shape[1], raxis), axis=0)
zlab = np.insert(get(ZLAB, nan_array(shape)), 0, np.full(shape[1], zaxis), axis=0)
jacobian = np.insert(
    get(JACOBIAN, nan_array(shape)), 0, np.full(shape[1], np.nan), axis=0
)

_g_norm = get(G_NORM, nan_array(shape[0]))
g_norm = np.insert(_g_norm, 0, _g_norm[0])
i_norm = np.insert(get(I_NORM, nan_array(shape[0])), 0, 0)
b_norm = np.insert(get(B_NORM, nan_array(shape)), 0, np.full(shape[1], 1), axis=0)

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

psip_norm_coord = xr.Variable(
    data=psip_norm,
    dims=["psip_norm"],
    attrs=dict(description="Poloidal flux coordinate", units="Normalized"),
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
    data=b,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="[T]"),
)

jacobian_var = xr.Variable(
    data=jacobian,
    dims=["psip_norm", "theta"],
    attrs=dict(description="VMEC to Boozer Jacobian", units="[m/T]"),
)

rlab_var = xr.Variable(
    data=rlab,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Lab horizontal coordinate", units="[m]"),
)

zlab_var = xr.Variable(
    data=zlab,
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
    data=b_norm,
    dims=["psip_norm", "theta"],
    attrs=dict(description="Magnetic field strength", units="Normalized"),
)


# ======================= Perturbations


# ========================================================


COORDS = {
    # "psip_norm": psip_norm_coord,
    "theta": theta_coord,
    "m": m_coord,
    "n": n_coord,
    "psi_norm": psi_norm_coord,
    "r_norm": r_norm_coord,
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
    "g_norm": g_norm_var,
    "i_norm": i_norm_var,
    "b_norm": b_norm_var,
}

description = "Converted from Falcon output"

ATTRS = {
    "description": description,
    "date": str(datetime.now().astimezone().replace(microsecond=0).isoformat()),
    "script": str(Path(__file__).stem),
    "version": str(VERSION),
}

dataset = xr.Dataset(
    data_vars=VARIABLES,
    coords=COORDS,
    attrs=ATTRS,
)

dataset.to_netcdf(OUTPUT)
GREEN = "\033[92m"
print(f"{GREEN}Stored dataset at '{OUTPUT.absolute()}'")

dataset.close()
raise SystemExit
