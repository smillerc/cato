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
init_pressure = np.float64(1e9) * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 3.5e13 * ureg("barye")
vacuum_density = 0.001 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

# Mesh
y_thickness = 1.0 * ureg("um")
print(f"y_thickness: {y_thickness.to('um')}")
dy = None  # will use smallest dx if None

layer_thicknesses = [.5, .5] * ureg("um")
layer_spacing = ["constant", "constant"]
layer_resolution = [20, 20] * ureg("1/um")

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)

layer_density = [shell_density, vacuum_density]
layer_u = [0, -vacuum_u.m] * ureg("cm/s")
#layer_u = [0, 0] * ureg("cm/s")
layer_v = [0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure]

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
    spacing_scale_factor=1.0,
    input_file="input.ini",
)

# Save to file
write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
# fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

# vc = ax1.pcolormesh(
#     domain["x"].to("um").m,
#     domain["y"].to("um").m,
#     domain["rho"].m,
#     edgecolor="k",
#     lw=0.001,
#     cmap="RdBu",
#     antialiased=True,
# )
# fig.colorbar(vc, ax=ax1, label="Density")
# ax1.set_xlabel("X")
# ax1.set_ylabel("Y")
# ax1.axis("equal")
# plt.show()
