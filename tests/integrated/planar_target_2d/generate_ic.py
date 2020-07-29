# -*- coding: utf-8 -*-
"""Make a simple 2D planar layered target for ICF-like implosions
and apply a perturbation at the ice-ablator interface"""
import matplotlib.pyplot as plt
from configparser import ConfigParser
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import *

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

# Physics
gamma = 5.0 / 3.0
init_pressure = np.float64(1e9) * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.005 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

# Mesh
vacuum_feathering = 1.1
y_thickness = 0.2 * ureg("um")
dy = None  # will use smallest dx if None

interface_loc = 20.0
layer_thicknesses = [interface_loc, 10, 3] * ureg("um")
layer_spacing = ["constant", "constant", "constant"]
layer_resolution = [20, 20, 20] * ureg("1/um")

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)

layer_density = [ice_density, shell_density, vacuum_density]
layer_u = [0, 0, 0] * ureg("cm/s")
layer_v = [0, 0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure, vacuum_pressure]

domain = make_2d_layered_grid(
    layer_thicknesses,
    layer_n_cells,
    layer_density,
    layer_u,
    layer_v,
    layer_pressure,
    y_thickness,
    dy=dy,
    layer_spacing=layer_spacing,
    spacing_scale_factor=vacuum_feathering,
    n_ghost_layers=2,
)

x = domain["xc"].to("cm").m
y = domain["yc"].to("cm").m

# Perturbation
x0 = ((interface_loc) * ureg("um")).to("cm").m
y0 = (y_thickness / 2.0).to("cm").m
# y0 = 0
fwhm = (1.5 * ureg("um")).to("cm").m
gaussian_order = 1
perturbation_frac = -0.5

# 2D
perturbation = (-perturbation_frac) * np.exp(
    -(((((x - x0) ** 2) / fwhm ** 2) + (((y - y0) ** 2) / fwhm ** 2)) ** gaussian_order)
).astype(np.float64)

# Cutoff the perturbation below 1e-6 otherwise there
# will be weird noise issues
if perturbation_frac < 0:
    perturbation[perturbation < 1e-3] = 0.0
else:
    perturbation[perturbation > 1e-3] = 0.0

perturbation = 1.0 - perturbation
# domain["rho"] = perturbation * domain["rho"]

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].to("um").m,
    domain["y"].to("um").m,
    domain["rho"].m,
    edgecolor="k",
    lw=0.001,
    cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
plt.show()
