# -*- coding: utf-8 -*-
"""
A simple script to view the results from the simulation
"""
import matplotlib.pyplot as plt
import h5py
import argparse

parser = argparse.ArgumentParser(description="View simulation results")
parser.add_argument("filename", type=str, help="Filename of the .hdf5 file to open")

args = parser.parse_args()

index_to_name = ["density", "u", "v", "pressure"]
data = {}

var_list = ["density", "x_velocity", "y_velocity", "pressure"]

state_list = [
    "leftright_midpoints_reference_state",
    "downup_midpoints_reference_state",
    "corner_reference_state",
    "evolved_leftright_midpoints_state",
    "evolved_downup_midpoints_state",
    "evolved_corner_state",
]
recon_state = ["reconstructed_state"]

with h5py.File(args.filename, "r") as h5:
    # Transpose to match the index convention within the code
    for var in state_list + var_list + recon_state:
        try:
            data[var] = h5[f"/{var}"][()].T
        except Exception:
            pass

# Plot the primitive quantities
fix, axes = plt.subplots(1, 4, figsize=(20, 8))
for i, name in enumerate(var_list):
    p = axes[i].imshow(data[name], origin="lower")
    plt.colorbar(p, ax=axes[i], fraction=0.046, pad=0.04)
    axes[i].set_title(name)
plt.suptitle("Primitive Quantities")
plt.tight_layout()
plt.show()

# Plot the reference and evolved state
for name in state_list:
    values = data[name]
    fix, axes = plt.subplots(1, 4, figsize=(20, 8))
    for i in range(4):
        p = axes[i].imshow(values[i, :, :], origin="lower")
        plt.colorbar(p, ax=axes[i], fraction=0.046, pad=0.04)
        axes[i].set_title(index_to_name[i])
    plt.suptitle(name)
    plt.tight_layout()
    plt.show()
