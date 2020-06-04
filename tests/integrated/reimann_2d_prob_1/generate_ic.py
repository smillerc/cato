# -*- coding: utf-8 -*-
"""Make a 2d reimann test problem"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
domain = make_uniform_grid(n_cells=(200, 200), xrange=(-1, 1), yrange=(-1, 1))

x = domain["xc"].m
y = domain["yc"].m
rho = domain["rho"].m
u = domain["u"].m
v = domain["v"].m
p = domain["p"].m

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.0
domain["v"] = domain["v"] * 0.0

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if x[i, j] >= 0 and y[i, j] >= 0:
            rho[i, j] = 1.1
            u[i, j] = 0.0
            v[i, j] = 0.0
            p[i, j] = 1.1
        elif x[i, j] >= 0 and y[i, j] <= 0:
            rho[i, j] = 0.5065
            u[i, j] = 0.0
            v[i, j] = 0.8939
            p[i, j] = 0.35
        elif x[i, j] <= 0 and y[i, j] >= 0:
            rho[i, j] = 0.5065
            u[i, j] = 0.8939
            v[i, j] = 0.0
            p[i, j] = 0.35
        elif x[i, j] <= 0 and y[i, j] <= 0:
            rho[i, j] = 1.1
            u[i, j] = 0.8939
            v[i, j] = 0.8939
            p[i, j] = 1.1

domain["rho"] = rho * ureg("g/cc")
domain["u"] = u * ureg("cm/s")
domain["v"] = v * ureg("cm/s")
domain["p"] = p * ureg("barye")

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    domain["p"].m,
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Pressure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
