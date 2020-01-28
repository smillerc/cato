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
domain = make_uniform_grid(n_cells=(500, 500), xrange=(-10, 10), yrange=(-10, 10))

E = 1
dr = 1
gamma = 1.4


# Set the initial conditions
def init_pressure(E, dr):
    nu = 1  # planar
    return (3 * (gamma - 1) * E) / ((nu + 1) * np.pi * dr ** nu)


rho = domain["rho"]
x = domain["xc"]
y = domain["yc"]

# Make pressure a centered gaussian with surrounding pressure of 1.0
domain["p"] = 100 * np.exp(-(x ** 2 + y ** 2)) + 1

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if (x[i, j] ** 2 + y[i, j] ** 2) > dr ** 2:
            domain["p"][i, j] = 0.1
        else:
            domain["p"][i, j] = 10  # init_pressure(E, dr)

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
