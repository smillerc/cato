# -*- coding: utf-8 -*-
"""Make the layered target test grid"""
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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

cells_per_micron = 20
gamma = 5.0 / 3.0
rho_shell = 1.0
vacuum_press = 1e8
vacuum_dens = 1.0e-5
vacuum_cs = np.sqrt(gamma * vacuum_press / vacuum_dens)
vacuum_mach = 0.0

vacuum_vel = vacuum_mach * vacuum_cs
vacuum_vel = np.sqrt(2.0 / (gamma + 1.0) * vacuum_press / rho_shell)

print(f"Vacuum Region Speed {vacuum_vel:.2e}")
layer_thicknesses = [1, 39, 10, 20] * ureg("um")
layer_n_cells = (layer_thicknesses.m * cells_per_micron).astype(np.int)
layer_n_cells[-1] = (layer_thicknesses[-1].m * 1).astype(np.int)
layer_density = [0.01, 0.25, 1.0, vacuum_dens] * ureg("g/cc")
layer_u = [0, 0, 0, vacuum_vel] * ureg("cm/s")
layer_v = [0, 0, 0, 0] * ureg("cm/s")
layer_pressure = [1e9, 1e9, 1e9, vacuum_press] * ureg("barye")

layered_target = make_1d_layered_grid(
    layer_thicknesses, layer_n_cells, layer_density, layer_u, layer_v, layer_pressure
)

last_shell_layer_index = sum(layer_thicknesses[:2]).m * cells_per_micron
print("last_shell_layer_index", last_shell_layer_index)
layered_target["rho"][last_shell_layer_index, :] *= 1.001

write_initial_hdf5(filename="layered_target", initial_condition_dict=layered_target)

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
for ax, v in zip([ax1, ax2], ["rho", "u"]):
    # vc = ax.plot(layered_target["xc"][:, 1], layered_target[v][:, 1])
    vc = ax.plot(layered_target[v][:, 1])
    ax.set_ylabel(v)
    ax.set_xlabel("X")

plt.show()
