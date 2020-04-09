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
vacuum_press = 5e8
vacuum_dens = 5.0e-4
vacuum_cs = np.sqrt(gamma * vacuum_press / vacuum_dens)
vacuum_mach = 0.0

vacuum_vel = vacuum_mach * vacuum_cs

print(f"Vacuum Region Speed {vacuum_vel:.2e}")
layer_thicknesses = [10, 5, 10] * ureg("um")
layer_n_cells = (layer_thicknesses.m * cells_per_micron).astype(np.int)
layer_density = [0.25, 1.0, vacuum_dens] * ureg("g/cc")
layer_u = [0, 0, vacuum_vel] * ureg("cm/s")
layer_v = [0, 0, 0] * ureg("cm/s")
layer_pressure = [1e9, 1e9, vacuum_press] * ureg("barye")

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
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.plot(layered_target["rho"][1, :])
ax1.set_xlabel("Array Index (i)")
ax1.set_ylabel("Density [g/cc]")

ax2 = ax1.twinx()
ax2.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
vc = ax2.plot(layered_target["u"][1, :], "--", color="k")
ax2.set_ylabel("Velocity [cm/s]")
plt.title("Initial Conditions")
plt.show()
