# -*- coding: utf-8 -*-
"""Make a 2d oblique contact discontinuity test problem"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
domain = make_uniform_grid(n_cells=(60, 50), xrange=(0, 0.6), yrange=(0, 0.5))

x = domain["xc"].m
y = domain["yc"].m
rho = domain["rho"].m

rho_L = 2.0
rho_R = 1.0

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.1
domain["v"] = domain["v"] * 0.1
domain["p"] = domain["p"] * 0.714

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if x[i, j] > y[i, j]:
            rho[i, j] = rho_L
        else:
            rho[i, j] = rho_R

domain["rho"] = rho * ureg("g/cc")

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    domain["rho"].m,
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
