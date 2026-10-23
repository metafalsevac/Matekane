#!/usr/bin/env python3
"""Plot Matekane output (optional helper; needs numpy + matplotlib).

    python3 scripts/plot.py output.txt      # 2D: altitude vs distance for each CL, with terrain
    python3 scripts/plot.py output3d.csv    # 3D: trajectory, attitude, alpha/beta/V, |q|-1

Figures are saved next to the input file as <name>_*.png (use --show to open windows).
"""
import os
import sys

import numpy as np
import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def ground(x):
    x = np.asarray(x, dtype=float)
    g = -7.54673e-12 * x**4 + 9.16044e-8 * x**3 - 1.8723e-4 * x**2 - 5.91546e-1 * x + 7.64963e3
    return np.where(x > 6233.58, 7480.58, g)


def plot2d(path, base):
    d = np.loadtxt(path, comments="#")
    fig, ax = plt.subplots(figsize=(10, 5))
    xs = np.linspace(0, 10000, 500)
    ax.fill_between(xs, 5800, ground(xs), color="tan", label="terrain")
    for cl in np.unique(d[:, 6]):
        r = d[d[:, 6] == cl]
        ax.plot(r[:, 0], r[:, 1], label=f"CL = {cl:g}")
    ax.set(xlabel="x (ft)", ylabel="altitude (ft)", title="Matekane takeoff (2D)", ylim=(5800, None))
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(base + "_2d.png", dpi=120, bbox_inches="tight")
    return [base + "_2d.png"]


def plot3d(path, base):
    d = np.genfromtxt(path, delimiter=",", names=True)
    out = []
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    xs = np.linspace(min(0, d["x_N"].min()), max(10000, d["x_N"].max()), 120)
    ys = np.linspace(min(d["y_E"].min(), -500), max(d["y_E"].max(), 500), 20)
    X, Y = np.meshgrid(xs, ys)
    ax.plot_surface(X, Y, ground(X), color="tan", alpha=0.4, linewidth=0)
    for cl in np.unique(d["CL_cmd"]):
        r = d[d["CL_cmd"] == cl]
        ax.plot(r["x_N"], r["y_E"], r["h"], label=f"CL = {cl:g}")
    ax.set(xlabel="North x (ft)", ylabel="East y (ft)", zlabel="altitude (ft)", title="Matekane takeoff (3D, NED)")
    ax.legend(fontsize=8)
    fig.savefig(base + "_3d_trajectory.png", dpi=120, bbox_inches="tight")
    out.append(base + "_3d_trajectory.png")

    fig, axs = plt.subplots(4, 1, figsize=(9, 10), sharex=True)
    for i, cl in enumerate(np.unique(d["CL_cmd"])):
        r = d[d["CL_cmd"] == cl]
        t, c = r["t"], f"C{i}"
        axs[0].plot(t, r["phi_deg"], color=c, label=f"roll CL={cl:g}")
        axs[0].plot(t, r["theta_deg"], "--", color=c, label=f"pitch CL={cl:g}")
        axs[1].plot(t, r["psi_deg"], color=c, label=f"CL={cl:g}")
        axs[2].plot(t, r["alpha_deg"], color=c, label=f"alpha CL={cl:g}")
        axs[2].plot(t, r["beta_deg"], "--", color=c, label=f"beta CL={cl:g}")
        axs[3].plot(t, r["qnorm_err"], color=c, label=f"CL={cl:g}")
    axs[0].set_ylabel("roll / pitch (deg)")
    axs[1].set_ylabel("yaw (deg)")
    axs[2].set_ylabel("alpha / beta (deg)")
    axs[3].set_ylabel("|q| - 1")
    axs[3].set_xlabel("t (s)")
    for a in axs:
        a.grid(True)
        a.legend(fontsize=7)
    fig.savefig(base + "_3d_attitude.png", dpi=120, bbox_inches="tight")
    out.append(base + "_3d_attitude.png")
    return out


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = args[0] if args else "output.txt"
    base = os.path.splitext(path)[0]
    with open(path) as f:
        first = f.readline()
    files = plot3d(path, base) if first.startswith("t,") else plot2d(path, base)
    print("wrote", ", ".join(files))
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
