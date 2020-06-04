# -*- coding: utf-8 -*-
"""Make the double periodic shear test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
double_shear = make_uniform_grid(
    n_cells=(200, 200), xrange=(-0.5, 0.5), yrange=(-0.5, 0.5)
)

# Set the initial conditions
double_shear["rho"] = double_shear["rho"] * 1.4
double_shear["p"] = double_shear["p"] * 4
u = double_shear["u"].m
v = double_shear["v"].m
x = double_shear["xc"].m
y = double_shear["yc"].m

# The u and v arrays depend on the location w/in the grid.
# Since they're cell-centered quantities, they need the location
# of the cell center (xc, yc)
v = np.sin(np.pi * (x + 1.5))
for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if y[i, j] < 0:
            u[i, j] = np.tanh(15 * (0.5 + y[i, j]))
        else:
            u[i, j] = np.tanh(15 * (0.5 - y[i, j]))

double_shear["u"] = u * ureg("cm/s")
double_shear["v"] = v * ureg("cm/s")

write_initial_hdf5(filename="double_shear", initial_condition_dict=double_shear)

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
