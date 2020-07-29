# -*- coding: utf-8 -*-
"""Make a simple 1D planar layered target for ICF-like implosions
and apply a perturbation at the ice-ablator interface"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse

sys.path.append(os.path.abspath("../../.."))
from pycato import *

# Physics
gamma = 5.0 / 3.0
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.01 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s").m

# Mesh
interface_loc = 20.0
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
    y_thickness=0,
    layer_spacing=layer_spacing,
    n_ghost_layers=2,
)

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)
