import os
import sys
import json
import matplotlib.pyplot as plt
import numpy as np
import pint
from pathlib import Path

sys.path.append(os.path.abspath("../../.."))
from pycato import *

two_pi = (2 * np.pi) * ureg("radian")

# Read grid specs from a json file
with open("grid_spec.json", "r") as f:
    grid = json.load(f)

wavelength = grid["perturbation_wavelength_um"] * ureg("um")
k = two_pi / wavelength

if grid["full_wavelength"]:
    y_thickness = two_pi / k
else:
    y_thickness = (two_pi / k) / 2.0

if "cell/wavelength" in grid["y_spacing"]:
    # set y-resolution based on wavelength
    print("Applying cells/wavelength")
    cells_per_wavelength = float(grid["y_spacing"].split(" ")[0])

    wavelength_fraction = y_thickness / wavelength
    cells_per_micron_y = int(cells_per_wavelength * wavelength_fraction.m)
    dy = y_thickness / cells_per_micron_y
else:
    print("Applying cells/micron")
    cells_per_micron_y = float(grid["y_spacing"].split(" ")[0]) * ureg("1 / um")
    dy = (1.0 / cells_per_micron_y.m) * ureg("um")

layer_thicknesses = grid["layer_thicknesses"] * ureg(grid["length_units"])
n_layers = len(grid["layer_thicknesses"])

if isinstance(grid["x_spacing"], list):
    raise Exception("Unable to handle a list for 'x_spacing' for now...")
else:
    layer_resolution = (
        [float(grid["x_spacing"].split(" ")[0])] * n_layers * ureg("1/um")
    )
    layer_spacing = ["constant"] * n_layers

layer_n_cells = np.round(
    (layer_thicknesses * layer_resolution).to_base_units()
).m.astype(int)


init_vac_vel = False
try:
    init_vac_vel = grid["init_vacuum_velocity"]
except:
    pass

gamma = grid["gamma"]
layer_density = grid["density"] * ureg(grid["density_units"])
layer_u = grid["x_velocity"] * ureg(grid["velocity_units"])
layer_v = grid["y_velocity"] * ureg(grid["velocity_units"])
layer_pressure = grid["pressure"] * ureg(grid["pressure_units"])

# Initialize the vacuum with a velocity
if init_vac_vel:
    print("Initializing a velocity in the vacuum layer")
    layer_u[-1] = -np.sqrt(
        2.0 / (gamma + 1.0) * layer_pressure[-1] / layer_density[-2]
    ).to("cm/s")
    print("Layer velocity", layer_u)

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
    spacing_scale_factor=None,
)

x = domain["xc"].to("um")
y = domain["yc"].to("um")

# Perturbation
x0 = grid["perturbation_x0"] * ureg("um")  # perturbation location
print(f"Perturbing at {x0}")

pert_x = np.exp(-k.m * ((x - x0).m) ** 2)
pert_x[pert_x < 1e-6] = 0.0

pert_y = -((1.0 - np.cos(k * y)) / 2.0).to_base_units().m
perturbation_loc = pert_x * pert_y

perturbation = np.abs(perturbation_loc * grid["perturbation_fraction"])
perturbation[perturbation < 1e-6] = 0.0

rho = domain["rho"]

# Cutoff the perturbation below 1e-6 otherwise there
# will be weird noise issues
rho = rho * (1.0 - perturbation)

domain["rho"] = rho

# Save to file
write_initial_hdf5(filename="initial_conditions", initial_condition_dict=domain)

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    domain["x"].to("um").m,
    domain["y"].to("um").m,
    domain["rho"].m,
    edgecolor="k",
    lw=0.1,
    cmap="RdBu_r",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.axis("equal")
# plt.show()
plt.savefig("ic.png")
