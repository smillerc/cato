# -*- coding: utf-8 -*-
"""Richtmyer-Meshkov Test Problem"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../../scripts"))
from generate_initial_grids import make_uniform_grid, write_initial_hdf5

# Make the empty grid
ni = 200
domain = make_uniform_grid(
    n_cells=(ni, ni * 5), xrange=(0.0, np.pi), yrange=(0.0, 4 * np.pi)
)

# Set the initial conditions
domain["p"][:, :] = 0.1
domain["u"][:, :] = 0.0
domain["v"][:, :] = 0.0
rho_heavy = 2.0
rho_light = 1.0

shock_pos = 11.0
interface = 3 * np.pi

for i in range(domain["xc"].shape[0]):
    for j in range(domain["xc"].shape[1]):
        x = domain["xc"][i, j]
        y = domain["yc"][i, j]

        if y > -0.5 * np.cos(2 * x) + interface:
            domain["rho"][i, j] = rho_heavy
        else:
            domain["rho"][i, j] = rho_light

        if y > shock_pos:
            domain["rho"][i, j] = 3.0
            domain["p"][i, j] = 1.0
            domain["v"][i, j] = -1.0

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="initial_conditions",
    initial_condition_dict=domain,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
fig, (ax1, ax2, ax3) = plt.subplots(figsize=(18, 8), nrows=1, ncols=3)

dens = ax1.pcolormesh(
    domain["x"], domain["y"], domain["rho"], edgecolor="k", lw=0.001, antialiased=True
)
fig.colorbar(dens, ax=ax1, label="Density")

uax = ax2.pcolormesh(
    domain["x"], domain["y"], domain["p"], edgecolor="k", lw=0.001, antialiased=True
)
fig.colorbar(uax, ax=ax2, label="Pressure")

vax = ax3.pcolormesh(
    domain["x"], domain["y"], domain["v"], edgecolor="k", lw=0.001, antialiased=True
)
fig.colorbar(vax, ax=ax3, label="Y Velocity")

for ax in [ax1, ax2, ax3]:
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.axis("equal")

plt.tight_layout()
plt.show()
