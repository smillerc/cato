# -*- coding: utf-8 -*-
"""Make the layered target test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pint

ureg = pint.UnitRegistry()

# sys.path.append(os.path.abspath('../../../'))
sys.path.append(os.path.abspath("../../../scripts"))
# from generate_initial_grids import make_uniform_grid, write_initial_hdf5
from generate_initial_grids import (
    make_1d_in_x_uniform_grid,
    write_initial_hdf5,
    make_1d_layered_grid,
)

cells_per_micron = 10

layer_thicknesses = [10, 40, 10, 30] * ureg("um")
layer_n_cells = (layer_thicknesses.m * cells_per_micron).astype(np.int)
# layer_n_cells = [10 , 50, 80, 10]
layer_density = [0.01, 0.24, 1.0, 1e-4, 1e-4] * ureg("g/cc")
layer_u = [0, 0, 0, 0, 0] * ureg("cm/s")
layer_v = [0, 0, 0, 0, 0] * ureg("cm/s")
layer_pressure = [1e9, 1e9, 1e9, 1e10, 1e10] * ureg("barye")

domain = make_1d_layered_grid(
    layer_thicknesses, layer_n_cells, layer_density, layer_u, layer_v, layer_pressure
)

bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

write_initial_hdf5(
    filename="initial_conditions",
    initial_condition_dict=domain,
    boundary_conditions_dict=bc_dict,
)

# Plot the results
# fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
# for ax, v in zip([ax1, ax2], ["rho", "p"]):
#     vc = ax.plot(domain["xc"][1, :], domain[v][1, :])
#     ax.set_ylabel(v)
#     ax.set_xlabel("X")

# plt.show()
