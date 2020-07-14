# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys

sys.path.append("../../..")
from pycato import *


def show_deltav(pert_ds, sym_ds, xrange, trange, title):
    """Find the difference in velocity between the datasets"""
    pert_v = pert_ds.x_velocity.squeeze().drop("y")
    sym_v = sym_ds.x_velocity.squeeze().drop("y")

    # Interpolate onto the other dataset
    delta_v = pert_v - sym_v.interp_like(pert_v)
    delta_v.name = "delta_v [kms/s]"

    xmin = xrange[0]
    xmax = xrange[1]
    ds_cutoff = delta_v.where((delta_v.x > xmin) & (delta_v.x < xmax), drop=True)

    plt.figure(figsize=(14, 6))
    ds_cutoff.max(axis=1).plot()
    plt.title(title)
    plt.savefig(f"{title}_deltav_hist.png")
    plt.show()

    plt.figure(figsize=(14, 6))

    for t in range(trange[0], trange[1] + trange[2], trange[2]):
        ds_cutoff[t].plot(label=f"t: {ds_cutoff[t].t.data:.2f} ns")
    plt.legend()
    plt.title(title)
    plt.savefig(f"{title}_slices.png")
    plt.show()

    return ds_cutoff


perturbed_fvleg = load_dataset(folder="FVLEG_perturbed")
symmetric_fvleg = load_dataset(folder="FVLEG_symmetric")

perturbed_ausm = load_dataset(folder="AUSM+-up_all_speed_perturbed")
symmetric_ausm = load_dataset(folder="AUSM+-up_all_speed_symmetric")

dv_ausm = show_deltav(
    pert_ds=perturbed_ausm,
    sym_ds=symmetric_ausm,
    xrange=[70, 80],
    trange=(600, 800, 50),
    title="AUSM+-up",
)
dv_fvleg = show_deltav(
    pert_ds=perturbed_fvleg,
    sym_ds=symmetric_fvleg,
    xrange=[70, 80],
    trange=(600, 800, 50),
    title="FVLEG",
)

plt.figure(figsize=(14, 8))
dv_ausm.max(axis=1).plot(label="AUSM+-up all-speed")
dv_fvleg.max(axis=1).plot(label="FVLEG")
plt.xlabel("Time [ns]")
plt.ylabel("Delta V (perturbed-symmetric) [km/s]")
plt.legend()
plt.savefig("growth_comparison.png")
plt.show()
