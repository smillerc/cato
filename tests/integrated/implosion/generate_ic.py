# -*- coding: utf-8 -*-
"""Make the implosion initial conditions"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import *

# Make the empty grid
domain = make_uniform_grid(
    input_file="input.ini",
    n_cells=(201, 201),
    xrange=(0.0, 0.3),
    yrange=(0.0, 0.3),
)

x = domain["xc"].m
y = domain["yc"].m
rho = domain["rho"].m
u = domain["u"].m
v = domain["v"].m
p = domain["p"].m

rho[:,:] = 1.0
u[:,:] = 0.0
v[:,:] = 0.0
p[:,:] = 1.0

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if x[i, j] + y[i, j] < 0.15:
            rho[i, j] = 0.125
            p[i, j] = 0.14


domain["rho"] = rho * ureg("g/cc")
domain["u"] = u * ureg("cm/s")
domain["v"] = v * ureg("cm/s")
domain["p"] = p * ureg("barye")

write_initial_hdf5(filename="ic", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].m,
    domain["y"].m,
    domain["p"].m,
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Pressure")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
