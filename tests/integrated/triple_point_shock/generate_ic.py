# -*- coding: utf-8 -*-
"""
Make the triple-point shock interaction initial conditions
See https://computing.llnl.gov/projects/blast/triple-point-shock-interaction for details
"""
import matplotlib.pyplot as plt
from configparser import ConfigParser
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import *

# Read the input file and make sure the spatial order is consistent
config = ConfigParser()
config.read("input.ini")
config.sections()
edge_interp = config["scheme"]["limiter"]
edge_interp = edge_interp.strip("'").strip('"')

if edge_interp in ["TVD5", "MLP5"]:
    n_ghost_layers = 3
else:
    n_ghost_layers = 2

# Make the empty grid
domain = make_uniform_grid(
    n_cells=(500, 300), xrange=(0, 7), yrange=(0, 3), n_ghost_layers=n_ghost_layers,
)

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
        if x[i, j] <= 1.0:
            rho[i, j] = 1.0
            p[i, j] = 1.0
        elif y[i, j] > 1.5:
            rho[i, j] = 0.125
            p[i, j] = 0.1
        else:
            rho[i, j] = 0.1
            p[i, j] = 0.1

domain["rho"] = rho * ureg("g/cc")
domain["p"] = p * ureg("barye")

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
fig.colorbar(vc, ax=ax1, label="Pressure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
