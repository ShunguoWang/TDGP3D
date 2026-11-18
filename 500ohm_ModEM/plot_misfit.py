#!/usr/bin/env python3
"""Plot data misfit of a ModEM inversion.

Usage: plot_misfit.py observed.dat predicted.dat [start_predicted.dat] [--out prefix]

Writes
  <out>_misfit_maps.png     normalised RMS per site, per data type and period
  <out>_misfit_summary.png  RMS per component and period
  <out>_profile.png         observed vs predicted along the y = 0 profile
                            (and the starting-model response, if given)
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

Z_COMPS = ["ZXX", "ZXY", "ZYX", "ZYY"]
T_COMPS = ["TX", "TY"]


def read_dat(path):
    """{(period, site, comp): (complex value, error, x, y)} from a ModEM data file."""
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) != 11 or line.lstrip().startswith(("#", ">")):
                continue
            try:
                per = float(p[0])
                x, y = float(p[4]), float(p[5])
                val = complex(float(p[8]), float(p[9]))
                err = float(p[10])
            except ValueError:
                continue
            out[(round(per, 8), p[1], p[7])] = (val, err, x, y)
    return out


def normalised_residuals(obs, pred):
    """Rows of (period, site, comp, x, y, r_re, r_im) with r = (obs - pred) / err."""
    rows = []
    for key, (v, e, x, y) in obs.items():
        if key not in pred:
            continue
        r = (v - pred[key][0]) / e
        rows.append((key[0], key[1], key[2], x, y, r.real, r.imag))
    return rows


def rms(r):
    r = np.asarray(r)
    return float(np.sqrt(np.mean(r ** 2))) if r.size else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("obs")
    ap.add_argument("pred")
    ap.add_argument("start", nargs="?", help="predicted data of the starting model")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.splitext(a.pred)[0]
    name = os.path.basename(a.pred)

    obs, pred = read_dat(a.obs), read_dat(a.pred)
    start = read_dat(a.start) if a.start else None
    rows = normalised_residuals(obs, pred)
    periods = sorted({r[0] for r in rows})
    total = rms([v for r in rows for v in r[5:]])
    print(f"overall RMS {total:.3f} from {2 * len(rows)} real data")
    if start:
        srows = normalised_residuals(obs, start)
        print(f"starting-model RMS {rms([v for r in srows for v in r[5:]]):.3f}")

    # ---- maps: per-site RMS for Z and T at each period, plus all data
    groups = [("Impedance", Z_COMPS), ("Tipper", T_COMPS)]
    fig, axs = plt.subplots(2, len(periods) + 1, figsize=(4 * (len(periods) + 1), 7.6),
                            constrained_layout=True)
    vmax = 2.0
    for gi, (gname, comps) in enumerate(groups):
        for pi, per in enumerate(periods + [None]):
            sel = [r for r in rows if r[2] in comps and (per is None or r[0] == per)]
            sites = {}
            for r in sel:
                sites.setdefault(r[1], [r[3], r[4], []])[2] += [r[5], r[6]]
            xy = np.array([[s[0], s[1], rms(s[2])] for s in sites.values()])
            ax = axs[gi, pi]
            sc = ax.scatter(xy[:, 1] / 1e3, xy[:, 0] / 1e3, c=xy[:, 2], s=38, marker="s",
                            cmap="viridis", vmin=0, vmax=vmax)
            ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, ec="r", ls="--", lw=1))
            ax.set_aspect("equal")
            label = "all periods" if per is None else f"T = {per:g} s ({1 / per:g} Hz)"
            ax.set_title(f"{gname}, {label}\nRMS {rms([v for r in sel for v in r[5:]]):.2f}",
                         fontsize=10)
            ax.set_xlabel("y, east (km)")
            ax.set_ylabel("x, north (km)")
    cb = fig.colorbar(sc, ax=axs, shrink=0.8)
    cb.set_label("site RMS (normalised residual)")
    fig.suptitle(f"{name}: per-site misfit, overall RMS {total:.3f} "
                 f"(red circle: sphere outline at depth)")
    fig.savefig(out + "_misfit_maps.png", dpi=150)
    plt.close(fig)

    # ---- summary: RMS per component and period (final vs starting model)
    comps = Z_COMPS + T_COMPS
    fig, ax = plt.subplots(figsize=(9, 4.2), constrained_layout=True)
    width = 0.8 / len(periods)
    for pi, per in enumerate(periods):
        vals = [rms([v for r in rows if r[0] == per and r[2] == c for v in r[5:]])
                for c in comps]
        xs = np.arange(len(comps)) + (pi - (len(periods) - 1) / 2) * width
        ax.bar(xs, vals, width, label=f"{1 / per:g} Hz (final)")
        if start:
            svals = [rms([v for r in srows if r[0] == per and r[2] == c for v in r[5:]])
                     for c in comps]
            ax.plot(xs, svals, "k_", ms=14, mew=2,
                    label="starting model" if pi == 0 else None)
    ax.axhline(1, color="gray", ls="--", lw=1)
    ax.set_xticks(np.arange(len(comps)), comps)
    ax.set_ylabel("RMS")
    ax.set_title(f"{name}: RMS per component and frequency")
    ax.legend(ncol=2, fontsize=9)
    fig.savefig(out + "_misfit_summary.png", dpi=150)
    plt.close(fig)

    # ---- profile along y = 0: observed vs predicted
    show = [("ZXY", "real"), ("ZXY", "imag"), ("ZYX", "real"), ("ZYX", "imag"),
            ("TX", "real"), ("TX", "imag")]
    fig, axs = plt.subplots(len(show), len(periods), figsize=(4.2 * len(periods), 2.1 * len(show)),
                            sharex=True, constrained_layout=True)
    for pi, per in enumerate(periods):
        for si, (comp, part) in enumerate(show):
            keys = sorted((k for k, v in obs.items()
                           if k[0] == round(per, 8) and k[2] == comp and abs(v[3]) < 1.0),
                          key=lambda k: obs[k][2])
            x = np.array([obs[k][2] for k in keys]) / 1e3
            f = (lambda c: c.real) if part == "real" else (lambda c: c.imag)
            o = np.array([f(obs[k][0]) for k in keys])
            e = np.array([obs[k][1] for k in keys])
            ax = axs[si, pi]
            ax.errorbar(x, o, e, fmt="o", ms=3, color="k", lw=0.8, label="observed")
            ax.plot(x, [f(pred[k][0]) for k in keys], "r-", lw=1.5, label="final model")
            if start:
                ax.plot(x, [f(start[k][0]) for k in keys], "b--", lw=1, label="starting model")
            ax.axvspan(-0.5, 0.5, color="orange", alpha=0.12)
            if pi == 0:
                ax.set_ylabel(f"{part[:2].capitalize()} {comp}")
            if si == 0:
                ax.set_title(f"T = {per:g} s ({1 / per:g} Hz)")
            if si == len(show) - 1:
                ax.set_xlabel("x, north (km) at y = 0")
    axs[0, 0].legend(fontsize=8)
    fig.suptitle(f"{name}: profile y = 0 (shaded: sphere extent). "
                 r"Z in [mV/km]/[nT], T dimensionless")
    fig.savefig(out + "_profile.png", dpi=130)
    plt.close(fig)

    for per in periods:
        z = rms([v for r in rows if r[0] == per and r[2] in Z_COMPS for v in r[5:]])
        t = rms([v for r in rows if r[0] == per and r[2] in T_COMPS for v in r[5:]])
        print(f"{1 / per:>5g} Hz: Z RMS {z:.3f}, T RMS {t:.3f}")
    print("wrote", ", ".join(out + s for s in ("_misfit_maps.png", "_misfit_summary.png",
                                               "_profile.png")))


if __name__ == "__main__":
    main()
