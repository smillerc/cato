# -*- coding: utf-8 -*-
"""Kelvin-Helmholtz test problem
Initial conditions taken from
"A Well Posed Kelvin-Helmholtz Instability Test Problem", C.P. McNally, et. al
See https://arxiv.org/pdf/1111.1764.pdf or
https://iopscience.iop.org/article/10.1088/0067-0049/201/2/18/meta"""

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

from configparser import ConfigParser
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("../../.."))
from pycato import make_uniform_grid, write_initial_hdf5, ureg

# Read the input file and make sure the spatial order is consistent
config = ConfigParser()
config.read("input.ini")
config.sections()
edge_interp = config["scheme"]["edge_interpolation_scheme"]
edge_interp = edge_interp.strip("'").strip('"')

if edge_interp in ["TVD3", "TVD5", "MLP3", "MLP5"]:
    n_ghost_layers = 2
else:
    n_ghost_layers = 1

# Make the empty grid
domain = make_uniform_grid(
    n_cells=(200, 200),
    xrange=(0.0, 1.0),
    yrange=(0.0, 1.0),
    n_ghost_layers=n_ghost_layers,
)

# Set the initial conditions
domain["p"] = domain["p"] * 2.5
rho = domain["rho"].m
u = domain["u"].m
v = domain["v"].m
rho_1 = 1.0
rho_2 = 2.0
L = 0.025
U_1 = 0.5
U_2 = -0.5
U_m = (U_1 - U_2) / 2.0
rho_m = (rho_1 - rho_2) / 2.0

for i in range(domain["xc"].shape[0]):
    for j in range(domain["xc"].shape[1]):
        x = domain["xc"][i, j].m
        y = domain["yc"][i, j].m

        # Perturbed velocity
        v[i, j] = 0.01 * np.sin(4.0 * np.pi * x)

        if 0.25 > y >= 0.0:
            rho[i, j] = rho_1 - rho_m * np.exp((y - 0.25) / L)
            u[i, j] = U_1 - U_m * np.exp((y - 0.25) / L)
        elif 0.5 > y >= 0.25:
            rho[i, j] = rho_2 + rho_m * np.exp((-y + 0.25) / L)
            u[i, j] = U_2 + U_m * np.exp((-y + 0.25) / L)
        elif 0.75 > y >= 0.5:
            rho[i, j] = rho_2 + rho_m * np.exp(-(0.75 - y) / L)
            u[i, j] = U_2 + U_m * np.exp(-(0.75 - y) / L)
        elif 1.0 > y >= 0.75:
            rho[i, j] = rho_1 - rho_m * np.exp(-(y - 0.75) / L)
            u[i, j] = U_1 - U_m * np.exp(-(y - 0.75) / L)


write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
try:
    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(18, 8), nrows=1, ncols=3)

    dens = ax1.pcolormesh(
        domain["x"],
        domain["y"],
        domain["rho"],
        edgecolor="k",
        lw=0.001,
        antialiased=True,
    )
    fig.colorbar(dens, ax=ax1, label="Density")

    uax = ax2.pcolormesh(
        domain["x"], domain["y"], domain["u"], edgecolor="k", lw=0.001, antialiased=True
    )
    fig.colorbar(uax, ax=ax2, label="X Velocity")

    vax = ax3.pcolormesh(
        domain["x"], domain["y"], domain["v"], edgecolor="k", lw=0.001, antialiased=True
    )
    fig.colorbar(vax, ax=ax3, label="Y Velocity")

    for ax in [ax1, ax2, ax3]:
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis("equal")

    plt.tight_layout()
    plt.show()
except Exception:
    pass
