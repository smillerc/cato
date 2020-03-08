# -*- coding: utf-8 -*-
"""Make the Woodward-Collela Blastwave Test Problem"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../../scripts"))
from generate_initial_grids import make_1d_in_x_uniform_grid, write_initial_hdf5

# Make the empty grid
shock_tube = make_1d_in_x_uniform_grid(n_cells=1300, limits=(0, 1.0))

# Set the initial conditions
shock_tube["rho"][:, :] = 1.0
shock_tube["u"][:, :] = 0.0
shock_tube["v"][:, :] = 0.0
gamma = 1.4
x = shock_tube["xc"]
y = shock_tube["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        # Right State
        if x[i, j] >= 0.9:
            shock_tube["p"][i, j] = 100.0
        # Left State
        elif x[i, j] <= 0.1:
            shock_tube["p"][i, j] = 1000.0
        # Middle state
        else:
            shock_tube["p"][i, j] = 0.01

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="initial_conditions",
    initial_condition_dict=shock_tube,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
fig, ax = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)
pax = ax.twinx()
dens = ax.plot(
    shock_tube["xc"][1, :], shock_tube["rho"][1, :], color="k", label="Density"
)
press = pax.plot(
    shock_tube["xc"][1, :], shock_tube["p"][1, :], color="b", label="Pressure"
)

lns = dens + press
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.set_ylabel("Density [g/cc]")
pax.set_ylabel("Pressure [barye]")
ax.set_xlabel("X")

plt.show()
