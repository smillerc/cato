# -*- coding: utf-8 -*-
"""Make a simple 1D planar layered target for ICF-like implosions
and apply a perturbation at the ice-ablator interface"""
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
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.01 * ureg("g/cc")

# Mesh
interface_loc = 20.0
layer_thicknesses = [interface_loc, 10, 2] * ureg("um")
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
    y_thickness=0,
    layer_spacing=layer_spacing,
    n_ghost_layers=n_ghost_layers,
)

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

try:
    plt.figure(figsize=(10, 2))
    plt.pcolormesh(
        domain["x"].m,
        domain["y"].m,
        domain["rho"].m,
        ec="k",
        lw=0.1,
        antialiased=True,
        cmap="viridis",
    )
    plt.axis("equal")
    plt.show()
except Exception:
    pass
