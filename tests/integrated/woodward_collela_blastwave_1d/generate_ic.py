# -*- coding: utf-8 -*-
"""Make the Woodward-Collela Blastwave Test Problem"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_1d_in_x_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
ic = make_1d_in_x_uniform_grid(n_cells=500, limits=(0, 1.0))

# Set the initial conditions
gamma = 1.4
x = ic["xc"].m
y = ic["yc"].m
p = ic["p"].m

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        # Right State
        if x[i, j] >= 0.9:
            p[i, j] = 100.0
        # Left State
        elif x[i, j] <= 0.1:
            p[i, j] = 1000.0
        # Middle state
        else:
            p[i, j] = 0.01

ic["p"] = p * ureg("barye")

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=ic)

# Plot the results
fig, ax = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)
pax = ax.twinx()
dens = ax.plot(ic["xc"][:, 1].m, ic["rho"][:, 1].m, color="k", label="Density")
press = pax.plot(ic["xc"][:, 1].m, ic["p"][:, 1].m, color="b", label="Pressure")

lns = dens + press
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.set_ylabel("Density [g/cc]")
pax.set_ylabel("Pressure [barye]")
ax.set_xlabel("X")

plt.show()
