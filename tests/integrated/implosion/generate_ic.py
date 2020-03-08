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
domain = make_uniform_grid(n_cells=(300, 300), xrange=(-0.5, 0.5), yrange=(-0.5, 0.5))

# Set the initial conditions
domain["rho"] = domain["rho"] * 1.0
p0 = 1
domain["p"] = domain["p"] * p0
domain["x"] = domain["x"]
domain["y"] = domain["y"]
x = domain["xc"]
y = domain["yc"]

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.0
domain["v"] = domain["v"] * 0.0

r0 = 0.2
for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        radius = np.sqrt(x[i, j] ** 2 + y[i, j] ** 2)

        if r0 + 0.05 < radius < r0 + 0.1:
            domain["p"][i, j] = p0
            domain["rho"][i, j] = 0.25
        elif r0 + 0.1 < radius < r0 + 0.15:
            domain["p"][i, j] = p0
            domain["rho"][i, j] = 1.0
        elif r0 + 0.16 < radius < r0 + 0.18:
            domain["p"][i, j] = 10
            domain["rho"][i, j] = 0.1
        else:
            domain["p"][i, j] = p0
            domain["rho"][i, j] = 0.1

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="initial_conditions",
    initial_condition_dict=domain,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=1, ncols=2)

vc = ax1.pcolormesh(
    domain["x"],
    domain["y"],
    domain["rho"],
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")

uc = ax2.pcolormesh(
    domain["x"],
    domain["y"],
    domain["p"],
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
fig.colorbar(uc, ax=ax2, label="Pressure")

ax1.axis("equal")
ax2.axis("equal")
plt.show()
