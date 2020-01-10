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
layered_target = make_uniform_grid(n_nodes=(200, 2), xrange=(0, 1), yrange=(0, 0.1))

# Set the initial conditions
p0 = 1e-3
# layered_target["rho"] = layered_target["rho"] * 1.0
layered_target["p"] = layered_target["p"] * p0
# Zero velocity everywhere
layered_target["u"] = layered_target["u"] * 0.0
layered_target["v"] = layered_target["v"] * 0.0

x = layered_target["xc"]
y = layered_target["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if 0 < x[i, j] <= 0.1:
            layered_target["rho"][i, j] = 0.025
        elif 0.1 < x[i, j] <= 0.4:
            layered_target["rho"][i, j] = 0.5
        elif 0.4 < x[i, j] <= 0.6:
            layered_target["rho"][i, j] = 1.0
        else:
            layered_target["rho"][i, j] = 1e-3


bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="layered_target",
    initial_condition_dict=layered_target,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    layered_target["x"],
    layered_target["y"],
    layered_target["rho"],
    edgecolor="k",
    lw=0.1,
    # cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
# ax1.axis("equal")
plt.show()
