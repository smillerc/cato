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
double_shear = make_uniform_grid(
    n_nodes=(256, 256), xrange=(0, 2 * np.pi), yrange=(0, 2 * np.pi)
)

# Set the initial conditions
rho_0 = np.pi / 15
delta = 1
double_shear["rho"] = double_shear["rho"] * rho_0
double_shear["p"] = double_shear["p"] * 4.0

x = double_shear["xc"]
y = double_shear["yc"]

# The u and v arrays depend on the location w/in the grid.
# Since they're cell-centered quantities, they need the location
# of the cell center (xc, yc)
double_shear["v"] = delta * np.sin(x)
for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if y[i, j] <= np.pi:
            double_shear["u"][i, j] = np.tanh((y[i, j] - np.pi / 2) / rho_0)
        else:
            double_shear["u"][i, j] = np.tanh((1.5 * np.pi - y[i, j]) / rho_0)

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="double_shear",
    initial_condition_dict=double_shear,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=1, ncols=2)

vc = ax1.pcolormesh(
    double_shear["x"],
    double_shear["y"],
    double_shear["v"],
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Y Velocity")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")

uc = ax2.pcolormesh(
    double_shear["x"],
    double_shear["y"],
    double_shear["u"],
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
ax2.set_xlabel("X")
ax2.set_ylabel("Y")

fig.colorbar(uc, ax=ax2, label="X Velocity")
ax1.axis("equal")
ax2.axis("equal")
plt.show()
