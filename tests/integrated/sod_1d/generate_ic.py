# -*- coding: utf-8 -*-
"""Make the 1D Sod shocktube initial conditions"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

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
shock_tube = make_1d_in_x_uniform_grid(
    n_cells=500, limits=(0, 1.0), n_ghost_layers=n_ghost_layers
)

# Set the initial conditions
shock_tube["u"] = shock_tube["u"] * 0.0
shock_tube["v"] = shock_tube["v"] * 0.0
shock_tube["p"] = shock_tube["p"].m
shock_tube["rho"] = shock_tube["rho"].m
gamma = 1.4
x = shock_tube["xc"]
y = shock_tube["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        # Left State
        if x[i, j].m <= 0.5:
            shock_tube["p"][i, j] = 1.0
            shock_tube["rho"][i, j] = 1.0
        # Right state
        else:
            shock_tube["p"][i, j] = 0.1
            shock_tube["rho"][i, j] = 0.125

shock_tube["p"] = shock_tube["p"] * ureg("barye")
shock_tube["rho"] = shock_tube["rho"] * ureg("g/cc")
write_initial_hdf5(filename="shock_tube_1d", initial_condition_dict=shock_tube)

# Plot the results
try:
    # fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
    # for ax, v in zip([ax1, ax2], ["rho", "p"]):
    #     vc = ax.plot(shock_tube["xc"][:, 1], shock_tube[v][:, 1], "-o")
    #     ax.set_ylabel(v)
    #     ax.set_xlabel("X")
    # plt.show()

    plt.figure(figsize=(10, 2))
    plt.pcolormesh(
        shock_tube["x"].m,
        shock_tube["y"].m,
        shock_tube["rho"].m,
        ec="k",
        lw=0.1,
        antialiased=True,
        cmap="viridis",
    )
    plt.axis("equal")
    plt.show()
except Exception:
    pass
