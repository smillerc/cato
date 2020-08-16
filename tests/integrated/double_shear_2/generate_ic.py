# -*- coding: utf-8 -*-
"""Make the double periodic shear test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
domain = make_uniform_grid(
    n_cells=(256, 256), xrange=(0, 2 * np.pi), yrange=(0, 2 * np.pi)
)

# Set the initial conditions
rho_0 = np.pi / 15
delta = 1
domain["rho"] = domain["rho"] * rho_0
domain["p"] = domain["p"] * 4.0
x = domain["xc"].m
y = domain["yc"].m
u = domain["u"].m

# The u and v arrays depend on the location w/in the grid.
# Since they're cell-centered quantities, they need the location
# of the cell center (xc, yc)
v = delta * np.sin(x)
for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if y[i, j] <= np.pi:
            u[i, j] = np.tanh((y[i, j] - np.pi / 2) / rho_0)
        else:
            u[i, j] = np.tanh((1.5 * np.pi - y[i, j]) / rho_0)

domain["u"] = u * ureg("cm/s")
domain["v"] = v * ureg("cm/s")

write_initial_hdf5(filename="double_shear", initial_condition_dict=domain)

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=1, ncols=2)

vc = ax1.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    domain["v"].m,
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Y Velocity")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")

uc = ax2.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    domain["u"].m,
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
