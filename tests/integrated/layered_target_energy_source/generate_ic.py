# -*- coding: utf-8 -*-
"""Make the layered target test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pint

sys.path.append(os.path.abspath("../../../scripts"))
from generate_initial_grids import (
    make_1d_in_x_uniform_grid,
    write_initial_hdf5,
    make_1d_layered_grid,
)

ureg = pint.UnitRegistry()

cells_per_micron = 10

layer_thicknesses = [10, 40, 10, 30] * ureg("um")
layer_n_cells = (layer_thicknesses.m * cells_per_micron).astype(np.int)
layer_density = [0.01, 0.24, 1.0, 1e-6] * ureg("g/cc")
layer_u = [0, 0, 0, 0] * ureg("cm/s")
layer_v = [0, 0, 0, 0] * ureg("cm/s")
layer_pressure = [1e9, 1e9, 1e9, 1e8] * ureg("barye")

layered_target = make_1d_layered_grid(
    layer_thicknesses, layer_n_cells, layer_density, layer_u, layer_v, layer_pressure
)

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="layered_target",
    initial_condition_dict=layered_target,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
# fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

# vc = ax1.plot(layered_target["rho"][1, :], )
# plt.show()
