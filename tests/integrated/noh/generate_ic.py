# -*- coding: utf-8 -*-
"""
Generate initial conditions for the Noh Implosion Test
Reference: https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-17-28269
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Make the empty grid
domain = make_uniform_grid(n_cells=(50, 50), xrange=(0, 1), yrange=(0, 1))

x = domain["xc"].m
y = domain["yc"].m

# Radial velocity
V_r = -1.0
theta = np.arctan2(y, x)

rho = np.ones_like(domain["rho"])
p = np.ones_like(domain["rho"])
p = p * (2 / 3) * 1e-12
u = V_r * np.cos(theta)
v = V_r * np.sin(theta)

# Add units
domain["rho"] = rho * ureg("g/cc")
domain["p"] = p * ureg("barye")
domain["u"] = u * ureg("cm/s")
domain["v"] = v * ureg("cm/s")

write_initial_hdf5(filename="noh", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    np.sqrt(domain["u"].m ** 2 + domain["v"].m ** 2),
    edgecolor="k",
    lw=0.001,
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Radial Velocity")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
