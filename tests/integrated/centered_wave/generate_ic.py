# -*- coding: utf-8 -*-
"""Make the double periodic shear test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# sys.path.append(os.path.abspath('../../../'))
sys.path.append(os.path.abspath("../../../scripts"))
from generate_initial_grids import make_uniform_grid, write_initial_hdf5

# Make the empty grid
domain = make_uniform_grid(n_nodes=(120, 120), xrange=(-0.5, 0.5), yrange=(-0.5, 0.5))

# Set the initial conditions
domain["rho"] = domain["rho"] * 1.0
# domain["p"] = domain["p"] * .001
p0 = 1e-3
domain["x"] = domain["x"]
domain["y"] = domain["y"]
x = domain["xc"]
y = domain["yc"]

# Make pressure a centered gaussian with surrounding pressure of 1.0
# domain["p"] = np.exp(-(x**2 + y**2)) * 1.0e6 + 1e6# 1 atm
domain["p"] = 2 * np.exp(-((x ** 2) / 0.01 + (y ** 2) / 0.01)) + p0

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.0
domain["v"] = domain["v"] * 0.0

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="gaussian", initial_condition_dict=domain, boundary_conditions_dict=bc_dict
)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"],
    domain["y"],
    domain["p"],
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Pessure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
