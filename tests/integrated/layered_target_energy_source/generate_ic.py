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

# Physics
gamma = 5.0 / 3.0
init_pressure = 1e9 * ureg("barye")

void_density = 0.01 * ureg("g/cc")
void_thickness = 1 * ureg("um")

ice_density = 0.25 * ureg("g/cc")
ice_thickness = 39 * ureg("um")

shell_density = 1.0 * ureg("g/cc")
shell_thickness = 10 * ureg("um")

vacuum_pressure = 1e8 * ureg("barye")
vacuum_density = 1e-5 * ureg("g/cc")
vacuum_thickness = 20 * ureg("um")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

# Mesh
cells_per_micron = 15 * ureg("1/um")  # all but the corona
vacuum_feathering = 1.08
y_thickness = 15 * ureg("um")
# dy = 0.1 * ureg("um")
dy = None  # will use smallest dx if None

layer_thicknesses = [
    void_thickness.to("cm").m,
    ice_thickness.to("cm").m,
    shell_thickness.to("cm").m,
    vacuum_thickness.to("cm").m,
] * ureg("cm")

layer_spacing = ["constant", "constant", "constant", "linear"]
layer_n_cells = [int(l * cells_per_micron) for l in layer_thicknesses]
layer_density = [void_density, ice_density, shell_density, vacuum_density]
layer_u = [0, 0, 0, vacuum_u.to("cm/s").m] * ureg("cm/s")
layer_v = [0, 0, 0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure, init_pressure, vacuum_pressure]

layered_target = make_1d_layered_grid(
    layer_thicknesses,
    layer_n_cells,
    layer_density,
    layer_u,
    layer_v,
    layer_pressure,
    layer_spacing=layer_spacing,
    spacing_scale_factor=vacuum_feathering,
)

write_initial_hdf5(filename="layered_target", initial_condition_dict=layered_target)

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
for ax, v in zip([ax1, ax2], ["rho", "u"]):
    vc = ax.plot(layered_target[v][:, 1], "o-")
    ax.set_ylabel(v)
    ax.set_xlabel("Grid Index (i)")

plt.show()

# Plot the results
fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
for ax, v in zip([ax1, ax2], ["rho", "u"]):
    vc = ax.plot(layered_target["xc"][:, 1].to("um"), layered_target[v][:, 1], "o-")
    ax.set_ylabel(v)
    ax.set_xlabel("X")

plt.show()
