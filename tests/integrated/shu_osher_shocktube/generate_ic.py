# -*- coding: utf-8 -*-
"""Make the 1D Shu-Osher Shock Tube"""
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
epsilon = 0.2
shock_tube["v"] = shock_tube["v"] * 0.0
u = shock_tube["u"].m
p = shock_tube["p"].m
rho = shock_tube["rho"].m
gamma = 1.4
x = shock_tube["xc"].m
y = shock_tube["yc"].m

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        # Left State
        if x[i, j] < 1 / 8:
            u[i, j] = 2.629369
            p[i, j] = 10.3333
            rho[i, j] = 3.857143
        # Right state
        else:
            u[i, j] = 0.0
            p[i, j] = 1.0
            rho[i, j] = 1.0 + 0.2 * np.sin(8.0 * x[i, j] * 2.0 * np.pi)


write_initial_hdf5(filename="initial_conditions", initial_condition_dict=shock_tube)

# Plot the results
try:
    fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
    for ax, v in zip([ax1, ax2], ["rho", "p"]):
        vc = ax.plot(shock_tube["xc"][:, 1], shock_tube[v][:, 1], "-o")
        ax.set_ylabel(v)
        ax.set_xlabel("X")
    plt.show()
except Exception:
    pass
