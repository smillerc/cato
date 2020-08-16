# -*- coding: utf-8 -*-
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

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
edge_interp = config["scheme"]["limiter"]
edge_interp = edge_interp.strip("'").strip('"')

if edge_interp in ["TVD5", "MLP5"]:
    n_ghost_layers = 3
else:
    n_ghost_layers = 2

# Physics
gamma = 5.0 / 3.0
init_pressure = np.float64(1e9) * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.01 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

# Mesh
vacuum_feathering = 1.1
wavelength = 2 * ureg("um")
two_pi = (2 * np.pi) * ureg("radian")
k = two_pi / wavelength
print(f"k: {k}")
y_thickness = (two_pi / k) / 2.0
print(f"y_thickness: {y_thickness.to('um')}")
dy = None  # will use smallest dx if None

interface_loc = 15.0
layer_thicknesses = [interface_loc, 10, 2] * ureg("um")
layer_spacing = ["constant", "constant", "constant"]
layer_resolution = [25, 25, 25] * ureg("1/um")

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
    n_ghost_layers=n_ghost_layers,
)

x = domain["xc"].to("um")
y = domain["yc"].to("um")

# Perturbation
x0 = (interface_loc + 5) * ureg("um")  # perturbation location

pert_x = np.exp(-k.m * ((x - x0).m) ** 2)
pert_x[pert_x < 1e-4] = 0.0

pert_y = -((1.0 - np.cos(k * y)) / 2.0).to_base_units().m
perturbation_loc = pert_x * pert_y
perturbation_frac = 0.5

perturbation = np.abs(perturbation_loc * perturbation_frac)
perturbation[perturbation < 1e-4] = 0.0

domain["rho"] = domain["rho"] * (1.0 - perturbation)
# Cutoff the perturbation below 1e-6 otherwise there
# will be weird noise issues

# Save to file
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
