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
domain = make_uniform_grid(n_cells=(401, 401), xrange=(-1, 1), yrange=(-1, 1))

domain["x"] = domain["x"]
domain["y"] = domain["y"]
x = domain["xc"]
y = domain["yc"]

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.0
domain["v"] = domain["v"] * 0.0

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if x[i, j] >= 0 and y[i, j] >= 0:
            domain["rho"][i, j] = 1.1
            domain["u"][i, j] = 0.0
            domain["v"][i, j] = 0.0
            domain["p"][i, j] = 1.1
        elif x[i, j] >= 0 and y[i, j] <= 0:
            domain["rho"][i, j] = 0.5065
            domain["u"][i, j] = 0.0
            domain["v"][i, j] = 0.8939
            domain["p"][i, j] = 0.35
        elif x[i, j] <= 0 and y[i, j] >= 0:
            domain["rho"][i, j] = 0.5065
            domain["u"][i, j] = 0.8939
            domain["v"][i, j] = 0.0
            domain["p"][i, j] = 0.35
        elif x[i, j] <= 0 and y[i, j] <= 0:
            domain["rho"][i, j] = 1.1
            domain["u"][i, j] = 0.8939
            domain["v"][i, j] = 0.8939
            domain["p"][i, j] = 1.1

bc_dict = {
    "+x": "zero_gradient",
    "+y": "zero_gradient",
    "-x": "zero_gradient",
    "-y": "zero_gradient",
}

write_initial_hdf5(
    filename="initial_conditions",
    initial_condition_dict=domain,
    boundary_conditions_dict=bc_dict,
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
fig.colorbar(vc, ax=ax1, label="Pressure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
