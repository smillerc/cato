# -*- coding: utf-8 -*-
"""Make the Sedov explosion initial conditions"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

from configparser import ConfigParser
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_1d_in_x_uniform_grid, write_initial_hdf5, ureg

# Read the input file and make sure the spatial order is consistent
config = ConfigParser()
config.read("input.ini")
config.sections()
edge_interp = config["scheme"]["edge_interpolation_scheme"]
edge_interp = edge_interp.strip("'").strip('"')

if edge_interp in ["TVD3", "TVD5", "MLP3", "MLP5"]:
    n_ghost_layers = 2
else:
    n_ghost_layers = 1

# Make the empty grid
domain = make_uniform_grid(
    n_cells=(100, 100),
    xrange=(-0.25, 0.25),
    yrange=(-0.25, 0.25),
    n_ghost_layers=n_ghost_layers,
)

# Set the initial conditions
domain["rho"] = domain["rho"] * 0.001
# domain["p"] = domain["p"] * .001
p0 = 1e-3
domain["x"] = domain["x"]
domain["y"] = domain["y"]
x = domain["xc"]
y = domain["yc"]

# Make pressure a centered gaussian with surrounding pressure of 1.0
fwhm = 0.001
p_max = 10.0
p = p_max * np.exp(-((x.m ** 2) / fwhm + (y.m ** 2) / fwhm)) + p0
p = (p - p.max()) + p_max  # Make the max actually 10, since it's slightly off
domain["p"] = p * ureg(str(domain["p"].units))

# Zero velocity everywhere
domain["u"] = domain["u"] * 0.0
domain["v"] = domain["v"] * 0.0

write_initial_hdf5(filename="sedov", initial_condition_dict=domain)

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
