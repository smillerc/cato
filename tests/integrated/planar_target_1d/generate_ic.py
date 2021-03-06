# -*- coding: utf-8 -*-
"""Make a simple 1D planar layered target for ICF-like implosions
and apply a perturbation at the ice-ablator interface"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import *

# Physics
gamma = 1.6666666667
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")

vacuum_pressure = 1.0e12 * ureg("barye")
vacuum_density = 0.01 * ureg("g/cc")
gas_density = 0.001 * ureg("g/cc")

v_shell = np.sqrt(2.0 / (gamma + 1.0) * vacuum_pressure / shell_density).to("cm/s").m

layer_thicknesses = [60, 35, 8, .5] * ureg("um")
n_layers = len(layer_thicknesses)


layer_spacing = ["constant"] * n_layers
r = 10
layer_resolution = [r, r, r, r] * ureg("1/um")

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)

layer_density = [gas_density, ice_density, shell_density, vacuum_density]
layer_u = [0, 0, 0, -v_shell] * ureg("cm/s")
layer_v = [0, 0, 0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure, init_pressure, vacuum_pressure]

domain = make_2d_layered_grid(
    layer_thicknesses,
    layer_n_cells,
    layer_density,
    layer_u,
    layer_v,
    layer_pressure,
    y_thickness=0,
    layer_spacing=layer_spacing,
)


# apply_perturbations = False
# if apply_perturbations:
#     print("Applying perturbations to the initial conditions")
#     x = domain["xc"].to("cm").m
#     y = domain["yc"].to("cm").m

#     # Perturbation
#     x0 = (68 * ureg("um")).to("cm").m
#     fwhm = (0.5 * ureg("um")).to("cm").m
#     gaussian_order = 1
#     perturbation_frac = -0.01

#     # 1D
#     perturbation = (-perturbation_frac) * np.exp(
#         -((((x - x0) ** 2) / fwhm ** 2) ** gaussian_order)
#     )

#     # Cutoff the perturbation below 1e-6 otherwise there
#     # will be weird noise issues
#     perturbation[perturbation < 1e-6] = 0.0

#     perturbation = 1.0 - perturbation
#     domain["rho"] = perturbation * domain["rho"]

write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# try:
#     plt.figure(figsize=(10, 2))
#     plt.pcolormesh(
#         domain["x"].m,
#         domain["y"].m,
#         domain["rho"].m,
#         ec="k",
#         lw=0.1,
#         antialiased=True,
#         cmap="viridis",
#     )
#     plt.axis("equal")
#     plt.show()
# except Exception:
#     pass
