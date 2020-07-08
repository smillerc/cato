# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_1d_layered_grid, ureg, write_initial_hdf5

# Physics
gamma = 5.0 / 3.0
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 7e12 * ureg("barye")
vacuum_density = 0.005 * ureg("g/cc")
vacuum_u = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s")

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

domain = make_1d_layered_grid(
    layer_thicknesses,
    layer_n_cells,
    layer_density,
    layer_u,
    layer_v,
    layer_pressure,
    layer_spacing=layer_spacing,
)

# Write the symmetric case
write_initial_hdf5(filename="symmetric", initial_condition_dict=domain)

x = domain["xc"].to("cm").m
y = domain["yc"].to("cm").m

# Perturbation
x0 = (interface_loc * ureg("um")).to("cm").m
fwhm = (0.5 * ureg("um")).to("cm").m
gaussian_order = 1
perturbation_frac = -0.01

# 1D
perturbation = (-perturbation_frac) * np.exp(
    -((((x - x0) ** 2) / fwhm ** 2) ** gaussian_order)
)

# Cutoff the perturbation below 1e-6 otherwise there
# will be weird noise issues
perturbation[perturbation < 1e-6] = 0.0

perturbation = 1.0 - perturbation
domain["rho"] = perturbation * domain["rho"]

write_initial_hdf5(filename="perturbed", initial_condition_dict=domain)

# Plot the results
# try:
# fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

# plt.plot(
#     domain["xc"][:, 1].to("um").m, domain["rho"][:, 1].m,
# )
# ax1.set_xlabel("X")
# ax1.set_ylabel("Density")
# plt.show()
# except Exception:
#     pass
