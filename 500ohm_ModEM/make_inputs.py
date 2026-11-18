#!/usr/bin/env python3
"""Build the ModEM 3D inversion input files for the Sphere data set.

Writes, next to this script:
  Sphere_prior.ws   starting/prior model (WS format, natural log resistivity)
  Sphere.cov        model covariance (smoothing) file
  control.inv       NLCG inversion control file
  control.fwd       forward solver control file

Mesh design (data: 361 sites, 200 m spacing, x,y = -1800..1800 m;
f = 10, 1, 0.1 Hz; 10 ohm-m sphere r = 500 m at 1 km depth in 1000 ohm-m):
  * 100 m core cells, sites at cell centres (2 cells per site spacing)
  * 15 padding cells per side, growth 1.5 -> ~130 km (> 2.5 skin depths at
    0.1 Hz in 1000 ohm-m, delta = 50 km)
  * 20 m top layer growing 1.1x to 50 m, held at 50 m down to 1.7 km
    (sphere spans 0.5-1.5 km), then 1.1x to ~250 m and 1.3x to > 150 km

Usage: make_inputs.py [rho]   uniform starting resistivity in ohm-m
                              (default: median apparent resistivity of data)
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data.dat")

CORE_DX = 100.0
CORE_HALF = 2150.0   # core spans -2150..2150 m -> 43 cells
NPAD = 15
PAD_FACTOR = 1.5

DZ_TOP = 20.0
DZ_FINE = 50.0       # constant dz through the target depth range
Z_FINE_BOTTOM = 1700.0
DZ_FACTOR_SHALLOW = 1.1
DZ_SHALLOW_MAX = 250.0
DZ_FACTOR_DEEP = 1.3
Z_BOTTOM = 150e3

SMOOTH = 0.3


def horizontal_cells():
    ncore = int(round(2 * CORE_HALF / CORE_DX))
    pad = [round(CORE_DX * PAD_FACTOR ** (k + 1)) for k in range(NPAD)]
    return pad[::-1] + [CORE_DX] * ncore + pad


def vertical_cells():
    dz, z, out = DZ_TOP, 0.0, []
    while z < Z_BOTTOM:
        out.append(round(dz))
        z += out[-1]
        if z < Z_FINE_BOTTOM:
            dz = min(dz * DZ_FACTOR_SHALLOW, DZ_FINE)
        else:
            dz *= DZ_FACTOR_SHALLOW if dz < DZ_SHALLOW_MAX else DZ_FACTOR_DEEP
    return out


def median_app_res():
    """Median Berdichevsky (det) apparent resistivity of the impedance data."""
    z = {}
    with open(DATA) as f:
        for line in f:
            p = line.split()
            if len(p) == 11 and p[7] in ("ZXY", "ZYX"):
                z.setdefault((p[0], p[1]), {})[p[7]] = complex(float(p[8]), float(p[9]))
    rho = []
    for (per, _), c in z.items():
        if len(c) == 2:
            # |Z| in [mV/km]/[nT]: rho_a = 0.2 T |Z|^2
            rho.append(0.2 * float(per) * abs(c["ZXY"] * c["ZYX"]))
    rho.sort()
    return rho[len(rho) // 2]


def fmt_rows(vals, per_line=10, fmt="{:g}"):
    return "\n".join(" ".join(fmt.format(v) for v in vals[i:i + per_line])
                     for i in range(0, len(vals), per_line))


def main():
    dx = horizontal_cells()
    dy = list(dx)
    dz = vertical_cells()
    nx, ny, nz = len(dx), len(dy), len(dz)

    if len(sys.argv) > 1:
        rho0 = float(sys.argv[1])  # starting resistivity given on the command line
    else:
        rho0 = float(f"{median_app_res():.1g}")  # one significant figure
    lnrho = math.log(rho0)

    with open(os.path.join(HERE, "Sphere_prior.ws"), "w") as f:
        f.write(f"# Sphere prior: {rho0:g} ohm-m half-space, {nx}x{ny}x{nz}\n")
        f.write(f"{nx} {ny} {nz} 0 LOGE\n")
        f.write(fmt_rows(dx) + "\n")
        f.write(fmt_rows(dy) + "\n")
        f.write(fmt_rows(dz) + "\n")
        row = fmt_rows([lnrho] * nx, fmt="{:.5E}")
        for _ in range(nz):
            f.write("\n")
            for _ in range(ny):
                f.write(row + "\n")
        # origin (m): grid centre at the surface, then rotation (deg)
        f.write(f"{-sum(dx) / 2:g} {-sum(dy) / 2:g} 0\n0\n")

    border = "+" + "-" * 77 + "+"
    header = [
        "This file defines model covariance for a recursive autoregression scheme.",
        "The model space may be divided into distinct areas using integer masks.",
        "Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between",
        "air, ocean and the rest of the model is turned off automatically. You can",
        "also define exceptions to override smoothing between any two model areas.",
        "To turn off smoothing set it to zero. This header is 16 lines long.",
        "1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)",
        "2. Smoothing in the X direction (NzEarth real values)",
        "3. Smoothing in the Y direction (NzEarth real values)",
        "4. Vertical smoothing (1 real value)",
        "5. Number of times the smoothing should be applied (1 integer >= 0)",
        "6. Number of exceptions (1 integer >= 0)",
        "7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 3 & 4)",
        "8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.",
    ]
    with open(os.path.join(HERE, "Sphere.cov"), "w") as f:
        f.write(border + "\n")
        for h in header:
            f.write(f"| {h:<76}|\n")
        f.write(border + "\n\n")
        f.write(f"{nx} {ny} {nz}\n\n")
        f.write(fmt_rows([SMOOTH] * nz) + "\n\n")
        f.write(fmt_rows([SMOOTH] * nz) + "\n\n")
        f.write(f"{SMOOTH}\n\n1\n\n0\n")

    with open(os.path.join(HERE, "control.inv"), "w") as f:
        f.write("Model and data output file name    : Sphere_NLCG\n"
                "Initial damping factor lambda      : 10.\n"
                "To update lambda divide by         : 10.\n"
                "Initial search step in model units : 1.\n"
                "Restart when rms diff is less than : 2.0e-3\n"
                "Exit search when rms is less than  : 1.05\n"
                "Exit when lambda is less than      : 1.0e-4\n"
                "Maximum number of iterations       : 200\n")

    with open(os.path.join(HERE, "control.fwd"), "w") as f:
        f.write("Number of QMR iters per divergence correction : 40\n"
                "Maximum number of divergence correction calls : 20\n"
                "Maximum number of divergence correction iters : 100\n"
                "Misfit tolerance for EM forward solver        : 1.0e-7\n"
                "Misfit tolerance for EM adjoint solver        : 1.0e-7\n"
                "Misfit tolerance for divergence correction    : 1.0e-5\n"
                # ModEM always reads line 7 as the nested-BC file; '#' = none
                "# Optional EM solution file name for nested BC: none\n"
                "Air layers mirror|fixed height|read from file : mirror\n"
                "Number of air layers and min top dz in km     : 10 3. 30.\n")

    print(f"grid {nx} x {ny} x {nz} = {nx * ny * nz} cells")
    print(f"x extent +/- {sum(dx) / 2 / 1e3:.1f} km, depth {sum(dz) / 1e3:.1f} km")
    print(f"dz: {dz}")
    print(f"median rho_a = {median_app_res():.0f} ohm-m -> prior {rho0:g} ohm-m")


if __name__ == "__main__":
    main()
