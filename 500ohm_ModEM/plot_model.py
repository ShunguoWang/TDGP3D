#!/usr/bin/env python3
"""Plot a vertical section (y = const) and a depth slice (z = const) of a
ModEM WS-format model, with the true sphere outlined for comparison.

Usage: plot_model.py model.rho [data.dat] [--y 0] [--z 993] [--out prefix]
"""
import argparse
import math
import os

import matplotlib.pyplot as plt
import numpy as np

# True model of the Sphere synthetic test
SPHERE_C = (0.0, 0.0, 1000.0)   # x (north), y (east), z (depth), m
SPHERE_R = 500.0


def read_ws(path):
    with open(path) as f:
        f.readline()
        head = f.readline().split()
        nx, ny, nz = map(int, head[:3])
        ptype = head[4] if len(head) > 4 else "LINEAR"
        vals = np.array(f.read().split(), dtype=float)
    dx, dy, dz = vals[:nx], vals[nx:nx + ny], vals[nx + ny:nx + ny + nz]
    n = nx * ny * nz
    v = vals[nx + ny + nz:nx + ny + nz + n]
    rest = vals[nx + ny + nz + n:]
    # file order: k slowest, then j, then i from Nx down to 1
    v = v.reshape(nz, ny, nx)[:, :, ::-1].transpose(2, 1, 0)  # -> (i, j, k)
    if ptype == "LOGE":
        rho = np.exp(v)
    elif ptype == "LOG10":
        rho = 10 ** v
    else:
        rho = v
    ox, oy = (rest[0], rest[1]) if len(rest) >= 2 else (-dx.sum() / 2, -dy.sum() / 2)
    xe = ox + np.concatenate([[0], np.cumsum(dx)])
    ye = oy + np.concatenate([[0], np.cumsum(dy)])
    ze = np.concatenate([[0], np.cumsum(dz)])
    return xe, ye, ze, rho


def read_sites(path):
    sites = set()
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 11 and not line.startswith(("#", ">")):
                sites.add((float(p[4]), float(p[5])))
    return np.array(sorted(sites))


def cell_index(edges, value):
    return int(np.clip(np.searchsorted(edges, value) - 1, 0, len(edges) - 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("data", nargs="?")
    ap.add_argument("--y", type=float, default=0.0, help="section position y (east), m")
    ap.add_argument("--z", type=float, default=993.0, help="slice depth, m")
    ap.add_argument("--xlim", type=float, default=2500.0, help="half-width shown, m")
    ap.add_argument("--zmax", type=float, default=3000.0, help="section depth shown, m")
    ap.add_argument("--clim", type=float, nargs=2, default=(1.0, 3.5), help="log10 rho range")
    ap.add_argument("--out", default=None, help="output prefix (default: model name)")
    a = ap.parse_args()

    xe, ye, ze, rho = read_ws(a.model)
    sites = read_sites(a.data) if a.data else None
    j = cell_index(ye, a.y)
    k = cell_index(ze, a.z)
    yc = 0.5 * (ye[j] + ye[j + 1])
    zc = 0.5 * (ze[k] + ze[k + 1])
    out = a.out or os.path.splitext(a.model)[0]
    name = os.path.basename(a.model)
    km = 1e-3
    norm = dict(vmin=a.clim[0], vmax=a.clim[1], cmap="jet_r", shading="flat")

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)

    # Vertical section x-z at y = const
    ax = axs[0]
    pc = ax.pcolormesh(xe * km, ze * km, np.log10(rho[:, j, :]).T, **norm)
    ax.set_xlim(-a.xlim * km, a.xlim * km)
    ax.set_ylim(a.zmax * km, 0)
    dyc = yc - SPHERE_C[1]
    if abs(dyc) < SPHERE_R:
        r = math.sqrt(SPHERE_R ** 2 - dyc ** 2)
        ax.add_patch(plt.Circle((SPHERE_C[0] * km, SPHERE_C[2] * km), r * km,
                                fill=False, ec="k", ls="--", lw=1.2))
    if sites is not None:
        s = sites[np.abs(sites[:, 1] - a.y) < 1.0]
        ax.plot(s[:, 0] * km, np.zeros(len(s)), "kv", ms=4, clip_on=False)
    ax.set_xlabel("x, north (km)")
    ax.set_ylabel("depth (km)")
    ax.set_title(f"Section y = {yc:.0f} m (cell {ye[j]:.0f} to {ye[j+1]:.0f} m)")
    ax.set_aspect("equal")

    # Depth slice x-y at z = const (x north up, y east right)
    ax = axs[1]
    ax.pcolormesh(ye * km, xe * km, np.log10(rho[:, :, k]), **norm)
    ax.set_xlim(-a.xlim * km, a.xlim * km)
    ax.set_ylim(-a.xlim * km, a.xlim * km)
    dzc = zc - SPHERE_C[2]
    if abs(dzc) < SPHERE_R:
        r = math.sqrt(SPHERE_R ** 2 - dzc ** 2)
        ax.add_patch(plt.Circle((SPHERE_C[1] * km, SPHERE_C[0] * km), r * km,
                                fill=False, ec="k", ls="--", lw=1.2))
    if sites is not None:
        ax.plot(sites[:, 1] * km, sites[:, 0] * km, "k.", ms=2)
    ax.axhline(0, color="w", lw=0.6, ls=":")
    ax.set_xlabel("y, east (km)")
    ax.set_ylabel("x, north (km)")
    ax.set_title(f"Depth slice z = {zc:.0f} m (layer {ze[k]:.0f} to {ze[k+1]:.0f} m)")
    ax.set_aspect("equal")

    cb = fig.colorbar(pc, ax=axs, shrink=0.85, pad=0.02)
    cb.set_label(r"log$_{10}$ resistivity ($\Omega$m)")
    fig.suptitle(f"{name}   (dashed: true 10 $\\Omega$m sphere; background 1000 $\\Omega$m)")
    fig.savefig(out + "_sections.png", dpi=150)

    imin = np.unravel_index(np.argmin(rho[:, :, :]), rho.shape)
    print(f"y section: cell j={j + 1}, centre y={yc:.1f} m")
    print(f"z slice: layer k={k + 1}, {ze[k]:.0f}-{ze[k + 1]:.0f} m, centre {zc:.1f} m")
    print(f"min rho {rho.min():.2f} ohm-m at x={0.5*(xe[imin[0]]+xe[imin[0]+1]):.0f}, "
          f"y={0.5*(ye[imin[1]]+ye[imin[1]+1]):.0f}, z={0.5*(ze[imin[2]]+ze[imin[2]+1]):.0f} m")
    print(f"slice rho range {rho[:, :, k].min():.1f}-{rho[:, :, k].max():.1f} ohm-m; "
          f"section range {rho[:, j, :].min():.1f}-{rho[:, j, :].max():.1f}")
    print("wrote", out + "_sections.png")


if __name__ == "__main__":
    main()
