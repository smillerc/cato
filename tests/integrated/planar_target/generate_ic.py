import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pint

sys.path.append(os.path.abspath("../../.."))
from pycato import (
    make_1d_in_x_uniform_grid,
    make_2d_layered_grid,
    make_1d_layered_grid,
    write_initial_hdf5,
)

ureg = pint.UnitRegistry()

# Physics
gamma = 5.0 / 3.0
init_pressure = np.float64(1e9) * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.005 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

# Mesh
vacuum_feathering = 1.1
y_thickness = .2 * ureg("um")
# dy = 0.1 * ureg("um")
dy = None  # will use smallest dx if None

interface_loc = 2.0
layer_thicknesses = [interface_loc, 10, 3] * ureg("um")
layer_spacing = ["constant", "constant", "constant"]
layer_resolution = [2, 20, 20] * ureg("1/um")

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
)

x = domain["xc"].to("cm").m
y = domain["yc"].to("cm").m

# Perturbation
x0 = (interface_loc * ureg("um")).to("cm").m
y0 = (y_thickness / 2.0).to("cm").m
fwhm = (0.5 * ureg("um")).to("cm").m
gaussian_order = 1
perturbation_frac = -0.01

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