# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "xarray",
#   "netcdf4",
#   "numpy",
#   "semver"
# ]
# ///

"""Converts a `npz` file to a `netcdf` file, with the expected
fields and variable names.
"""

import semver

VERSION = semver.Version.parse("0.1.0-alpha.1")

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "input",
    help="the input 'npz' file",
)
parser.add_argument(
    "output",
    help="the output 'netCDF' file",
)
parser.add_argument(
    "-q",
    "--quiet",
    help="do not print messages",
    action="store_true",
)
args = parser.parse_args()

# ==========================================================================

import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime
from warnings import warn


if args.quiet:
    import sys
    import io

    text_trap = io.StringIO()
    sys.stdout = text_trap

xr.set_options(
    display_max_rows=30,
    display_width=80,
    netcdf_engine_order=["netcdf4"],
)

INPUT = Path(args.input)
OUTPUT = Path(args.output)
npz = np.load(INPUT)

print(f"Loaded .npz file from '{INPUT.absolute()}'")

# NOTE: Some arrays seem to have endianness explicitly stored in the dtype
# ('<' or '>'). Convert to the native endianness ('=', default) and expected
# type, just in case.

# NOTE: We need to divide both fluxes by 2π.

baxis = npz.get("BB")[0, 0].astype("f8")  # [Tesla]
raxis = npz.get("Rmaj").astype("f8")  # [meters]
psip = npz.get("psipol").astype("f8") / (2 * np.pi)
psi = npz.get("psitor").astype("f8") / (2 * np.pi)
theta = npz.get("theta").astype("f8")
r = npz.get("r").astype("f8")
i = npz.get("I").astype("f8")
g = npz.get("g").astype("f8")
q = npz.get("q").astype("f8")
b = npz.get("BB").astype("f8")
rlab = npz.get("RR").astype("f8")
zlab = npz.get("ZZ").astype("f8")

# If no perturbations are found in the npz file, this will create empty
# Variables, which we can drop at the end.
default_coord = np.array([])
default_array = np.full((0, 0, len(psip)), np.nan)
m = npz.get("m", default_coord).astype("i8")
n = npz.get("n", default_coord).astype("i8")
alphas = npz.get("alphas", default_array).astype("f8") / (2 * np.pi) ** 2
phases = npz.get("phases", default_array).astype("f8")

print("Exctracted all variables from npz file")

# Normalize
psip_norm = psip / (baxis * raxis**2)
psi_norm = psi / (baxis * raxis**2)
r_norm = r / raxis
g_norm = g / (baxis * raxis)
i_norm = i / (baxis * raxis)
b_norm = b / baxis
alphas_norm = alphas / raxis


# ======================= Scalars


baxis_var = xr.Variable(
    data=baxis,
    dims=[],
    attrs=dict(description="The magnetic field strength on the axis `B0`", units="[T]"),
)

raxis_var = xr.Variable(
    data=raxis,
    dims=[],
    attrs=dict(description="The major radius `R`", units="[m]"),
)


# ======================= Coordinates

# Apart from `θ`, `ψp`, `ψ` and `r`, we also use the mode numbers `m` and
# `n` as coordinates. This way, the `a{m,n}(ψp)` 1D arrays can be easily
# accessed with `dataset.alphas[m, n]`

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

psi_var = xr.Variable(
    data=psi,
    dims=["psi_norm"],
    attrs=dict(description="Toroidal flux coordinate", units="[Tm]"),
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


# WARN: Replace `inf` and `nan` values that may appear with 0.
np.nan_to_num(alphas_var.to_numpy(), nan=0, posinf=0, neginf=0, copy=False)
np.nan_to_num(alphas_norm_var.to_numpy(), nan=0, posinf=0, neginf=0, copy=False)


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
    "psip": psip_var,
    "psi": psi_var,
    "r": r_var,
    "q": q_var,
    "g": g_var,
    "i": i_var,
    "b": b_var,
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

# Remove empty Variables in case of no perturbations
dataset = dataset.drop_vars(
    [
        dataset[item].name
        for item in dataset.coords | dataset.data_vars
        if dataset[item].size == 0
    ]
)

# Useful metadata to avoid confusion
dataset = dataset.assign_attrs(
    {
        "input file": str(INPUT),
        "output file": str(OUTPUT),
        "date": str(datetime.now().astimezone().replace(microsecond=0).isoformat()),
        "script": str(Path(__file__).stem),
        "convention version": str(VERSION),
    }
)

print("Created dataset")

# makes `ncdump -h` cleaner
for item in dataset.coords | dataset.data_vars:
    dataset[item].attrs["_FillValue"] = False


for item in dataset.coords | dataset.data_vars:
    if not np.all(np.isfinite(dataset[item].to_numpy())):
        warn(f"NaN or inf found in `{dataset[item].name}`")

dataset.close()
dataset.to_netcdf(OUTPUT)

print(f"Stored dataset at '{OUTPUT.absolute()}'")

raise SystemExit
