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


# Physics
gamma = 5.0 / 3.0
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")
vacuum_density = 1e-3 * ureg("g/cc")

# Mesh
y_thickness = 10 * ureg("um")
print(f"y_thickness: {y_thickness.to('um')}")
dy = None  # will use smallest dx if None

# Mesh
interface_loc = 10.0
layer_thicknesses = [interface_loc, 10, 2] * ureg("um")
layer_spacing = ["constant", "constant", "constant"]
res = 10
layer_resolution = [res, res, res] * ureg("1/um")

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)

layer_density = [ice_density, shell_density, vacuum_density]
layer_u = [0, 0, 0] * ureg("cm/s")
layer_v = [0, 0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure, init_pressure]

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
    spacing_scale_factor=None,
    input_file="input.ini",
)

x = domain["xc"].to("um")
y = domain["yc"].to("um")

# Perturbation
do_pert = False
if do_pert:
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
plot = False
if plot:
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
