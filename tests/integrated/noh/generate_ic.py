# -*- coding: utf-8 -*-
"""Generate initial conditions for the Noh Implosion Test"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# sys.path.append(os.path.abspath('../../../'))
sys.path.append(os.path.abspath("../../../scripts"))
from generate_initial_grids import make_uniform_grid, write_initial_hdf5

# Make the empty grid
domain = make_uniform_grid(n_cells=(200, 200), xrange=(0, 1), yrange=(0, 1))

x = domain["xc"]
y = domain["yc"]
V_r = -1.0  # [cm/s] inward radial velocity
theta = np.arctan2(y, x)

domain["rho"] = np.ones_like(domain["rho"])  # [g/cc] uniform density
domain["p"] = np.ones_like(domain["rho"])  # [barye] uniform density

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        radius = np.sqrt(x[i, j] ** 2 + y[i, j] ** 2)
        theta = np.arctan2(y[i, j], x[i, j])
        if radius <= 0.9:
            domain["u"][i, j] = V_r * np.cos(theta)
            domain["v"][i, j] = V_r * np.sin(theta)
        else:
            domain["u"][i, j] = 0.0
            domain["v"][i, j] = 0.0

bc_dict = {"+x": "periodic", "-x": "periodic", "+y": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="ic", initial_condition_dict=domain, boundary_conditions_dict=bc_dict
)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"],
    domain["y"],
    np.sqrt(domain["u"] ** 2 + domain["v"] ** 2),
    edgecolor="k",
    lw=0.001,
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Radial Velocity")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
