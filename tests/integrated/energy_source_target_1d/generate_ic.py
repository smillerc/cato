# -*- coding: utf-8 -*-
"""Make a simple 1D planar layered target for ICF-like implosions
and apply a perturbation at the ice-ablator interface"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse
from pathlib import Path

sys.path.append(os.path.abspath("../../.."))
from pycato import *

parser = argparse.ArgumentParser()
parser.add_argument(
    "perturbed",
    type=str,
    help="Write out the perturbed initial conditions",
    default="symmetric",
)
args = parser.parse_args()
apply_perturbations = False
if args.perturbed == "perturbed":
    apply_perturbations = True

# Physics
gamma = 5.0 / 3.0
fill_density = 1e-3 * ureg("g/cc")
init_pressure = 1e9 * ureg("barye")
ice_density = 0.25 * ureg("g/cc")
shell_density = 1.0 * ureg("g/cc")
vacuum_press = 10 * ureg('Mbar')
vacuum_density = 0.005 * ureg("g/cc")

# Mesh
layer_thicknesses = [30, 6, 3] * ureg("um")
interface_loc = np.cumsum(layer_thicknesses)[1]
gas_ice_interface_loc = layer_thicknesses[0]
layer_spacing = ["constant"] * len(layer_thicknesses)
cells_per_micron_x = 10
dy = None  # 1D

layer_resolution = [cells_per_micron_x, cells_per_micron_x, cells_per_micron_x,] * ureg(
    "1/um"
)

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)

v_shell = np.sqrt(2.0 / (gamma + 1.0) * vacuum_press / shell_density).to("cm/s").m

layer_density = [ice_density, shell_density, vacuum_density]
layer_u = [0, 0, 0] * ureg("cm/s")
layer_v = [0, 0, 0] * ureg("cm/s")
layer_pressure = [init_pressure, init_pressure, init_pressure]

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

if False:
    print("Applying perturbations to the initial conditions")
    x = domain["xc"].to("um")

    # Perturbation
    wavelength = 100 * ureg("um")
    two_pi = (2 * np.pi) * ureg("radian")
    k = two_pi / wavelength

    # perturbation location
    # x0 = interface_loc - 2 * ureg('um')
    x0 = gas_ice_interface_loc + 2 * ureg("um")

    pert_x = np.exp(-k.m * ((x - x0).m) ** 2)
    pert_x[pert_x < 1e-6] = 0.0
    perturbation_loc = pert_x
    perturbation_frac = 0.1

    perturbation = np.abs(perturbation_loc * perturbation_frac)
    perturbation[perturbation < 1e-6] = 0.0

    rho = domain["rho"]

    # Only apply the perturbation in the ice
    mask = np.where(rho.m <= 0.25)

    # Cutoff the perturbation below 1e-6 otherwise there
    # will be weird noise issues
    rho[mask] = rho[mask] * (1.0 - perturbation[mask])

    domain["rho"] = rho

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
