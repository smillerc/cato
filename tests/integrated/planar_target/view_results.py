# -*- coding: utf-8 -*-
"""
A simple script to view the results from the simulation
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append("../../..")
from pycato import load_2d_dataset

# # Load cato results
ds = load_2d_dataset("results", units="icf")
ds["mach_y"] = ds.y_velocity / ds.sound_speed
ds["mach_x"] = ds.x_velocity / ds.sound_speed

t = -1

max_mach_y = abs(ds.mach_y[t]).max().data
fig, (ax, ax1, ax2) = plt.subplots(nrows=3, sharex=True, ncols=1, figsize=(12, 8))

# Contour
ax.axis("equal")
ds.density[t].plot.pcolormesh(
    x="x", y="y", lw=0.1, antialiased=True, ax=ax, cmap="viridis"
)

ax1.axis("equal")
ds.mach_y[t].plot.pcolormesh(
    x="x", y="y", lw=0.1, antialiased=True, ax=ax1, cmap="RdBu_r"
)
ax1.set_title(f"Max Mach_y: {max_mach_y}")

rho_ax = ax2.twinx()
# 1D Lineout
for v in ["density", "x_velocity", "pressure"]:
    var = ds[v][t]
    var_dy = var[:, 2] - var[:, 1]
    var_dy.plot(x="x", ax=ax2, label=f"d/dy {v}")
ds.density[t][:, 1].plot(x="x", ax=rho_ax, label="Density", color="k")
# ds.pressure[-1][:,1].plot(x='x', ax=ax2, label='Pressure')
# ds.mach_x[-1][:,1].plot(x='x', ax=ax2, label='Mach_x')

ax2.legend()
ax2.set_ylabel("")
plt.tight_layout()
plt.savefig("results.png")

try:
    plt.show()
except Exception:
    pass
