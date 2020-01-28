# -*- coding: utf-8 -*-
"""Make the layered target test grid"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import pint

ureg = pint.UnitRegistry()

# sys.path.append(os.path.abspath('../../../'))
sys.path.append(os.path.abspath("../../../scripts"))
# from generate_initial_grids import make_uniform_grid, write_initial_hdf5
from generate_initial_grids import make_1d_in_x_uniform_grid, write_initial_hdf5

# Make the empty grid
layer_thicknesses = np.array([20, 30, 10, 40]) * ureg("um")
layer_densities = np.array([0.01, 0.25, 1.0, 1e-3]) * ureg("g/cc")
layer_pressures = np.array([1e6, 1e9, 1e9, 1e4]) * ureg("barye")
layered_target = make_1d_in_x_uniform_grid(
    n_cells=550, limits=(0, np.sum(layer_thicknesses).to("cm"))
)  # 100 microns

print(np.sum(layer_thicknesses).to("cm"))
# Set the initial conditions
p0 = 1e6  # baryes is about 1 atm
# layered_target["rho"] = layered_target["rho"] * 1.0
layered_target["p"] = layered_target["p"] * p0
# Zero velocity everywhere
layered_target["u"] = layered_target["u"] * 0.0
layered_target["v"] = layered_target["v"] * 0.0

x = layered_target["xc"]
y = layered_target["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        # print(x[i, j])
        if 0 < x[i, j] <= 0.001:  # inner
            # print('gas')
            layered_target["rho"][i, j] = layer_densities[0].m
            layered_target["p"][i, j] = layer_pressures[0].m

        elif 0.001 < x[i, j] <= 0.004:  # ice
            # print('ice')
            layered_target["rho"][i, j] = layer_densities[1].m
            layered_target["p"][i, j] = layer_pressures[1].m

        elif 0.004 < x[i, j] <= 0.007:  # shell
            # print('shell')
            layered_target["rho"][i, j] = layer_densities[2].m
            layered_target["p"][i, j] = layer_pressures[2].m

        else:  # vaccuum
            # print('vaccuum')
            layered_target["rho"][i, j] = layer_densities[3].m
            layered_target["p"][i, j] = layer_pressures[3].m

# layer = 0
# layer_thicknesses_cm = layer_thicknesses.to('cm').m
# for i in range(x.shape[1]):
#     r = x[0, i]
#     print(r, np.cumsum(layer_thicknesses[:layer + 1]))
# if r > np.cumsum(layer_thicknesses[:layer].m):
# layer += 1
# print(layer, r)
# if layer == 0:
#     layered_target["rho"][i, :] = layer_densities[0]
#     layered_target["p"][i, :] = layer_pressures[0]
# else:
#     layered_target["rho"][i, :] = layer_densities[layer]
#     layered_target["p"][i, :] = layer_pressures[layer]

# if 0 < r <= 0.001:  # inner
#     layered_target["rho"][i, :] = 0.01
#     layered_target["p"][i, :] = 1e6

# elif 0.001 < r <= 0.4:  # ice
#     layered_target["rho"][i, :] = 0.25
#     layered_target["p"][i, :] = 1e9

# elif 0.4 < r <= 0.7:  # shell
#     layered_target["rho"][i, :] = 1.0
#     layered_target["p"][i, :] = 1e9

# else:  # vaccuum
#     layered_target["rho"][i, :] = 1e-3
#     layered_target["p"][i, :] = 1e4

# bc_dict = {
#     "+x": "periodic",
#     "+y": "periodic",
#     "-x": "periodic",
#     "-y": "periodic"
# }

# write_initial_hdf5(
#     filename="layered_target",
#     initial_condition_dict=layered_target,
#     boundary_conditions_dict=bc_dict,
# )

# Plot the results
fig, (ax1) = plt.subplots(figsize=(18, 8), nrows=1, ncols=1)

vc = ax1.pcolormesh(
    layered_target["x"],
    layered_target["y"],
    layered_target["rho"],
    edgecolor="k",
    lw=0.1,
    # cmap="RdBu",
    antialiased=True,
)
fig.colorbar(vc, ax=ax1, label="Density")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
# ax1.axis("equal")
plt.show()
