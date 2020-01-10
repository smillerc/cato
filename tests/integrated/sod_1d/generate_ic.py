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
shock_tube = make_uniform_grid(n_nodes=(200, 2), xrange=(0, 1), yrange=(0, 1))

# Set the initial conditions
shock_tube["u"] = shock_tube["u"] * 0.0
shock_tube["v"] = shock_tube["v"] * 0.0

x = shock_tube["xc"]
y = shock_tube["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if x[i, j] <= 0.5:
            shock_tube["rho"][i, j] = 1.0
            shock_tube["p"][i, j] = 1.0
        else:
            shock_tube["rho"][i, j] = 0.125
            shock_tube["p"][i, j] = 0.1

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="shock_tube_1d",
    initial_condition_dict=shock_tube,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
# fig, (ax1, ax2, ax3) = plt.subplots(figsize=(18, 8), nrows=3, ncols=1)
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)

for ax, v in zip([ax1, ax2], ["rho", "p"]):
    vc = ax.pcolormesh(
        shock_tube["x"],
        shock_tube["y"],
        shock_tube[v],
        edgecolor="k",
        lw=0.1,
        antialiased=True,
    )
    fig.colorbar(vc, ax=ax, label=v)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

# ax3.plot(shock_tube["xc"][0,:], shock_tube['rho'][0,:], label='Density')
# ax3.plot(shock_tube["xc"][0,:], shock_tube['p'][0,:], label='Pressure')
# ax3.legend()
# plt.tight_layout()
plt.show()
