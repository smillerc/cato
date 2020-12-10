# -*- coding: utf-8 -*-
"""Make the 1D Sod shocktube initial conditions"""
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


xmax = 100
shell_width = 3.0
shock_tube = make_1d_in_x_uniform_grid(
    n_cells=5000, limits=(0, xmax), input_file="input.ini"
)

middle = xmax / 2.0
width = 0.1

# Set the initial conditions
shock_tube["u"] = shock_tube["u"] * 0.0
shock_tube["v"] = shock_tube["v"] * 0.0
shock_tube["p"] = shock_tube["p"].m
shock_tube["rho"] = shock_tube["rho"].m
gamma = 1.4
x = shock_tube["xc"]
y = shock_tube["yc"]

for i in range(y.shape[0]):
    for j in range(y.shape[1]):
        if  middle - shell_width <= x[i, j].m <= middle + shell_width:
            shock_tube["rho"][i, j] = 1.0
        else:
            shock_tube["rho"][i, j] = 0.25

        if middle - width <= x[i, j].m <= middle + width:
            shock_tube["p"][i, j] = 1e4
            shock_tube["rho"][i, j] = 0.1
        else:
            shock_tube["p"][i, j] = 0.1
# for i in range(y.shape[0]):
#     for j in range(y.shape[1]):
#         if  9.0 <= x[i, j].m <= 10.0:
#             shock_tube["rho"][i, j] = 1.0
#         else:
#             shock_tube["rho"][i, j] = 0.25

#         if x[i, j].m > 10.0:
#             shock_tube["p"][i, j] = 1e3
#             shock_tube["rho"][i, j] = 0.1
#         else:
#             shock_tube["p"][i, j] = 0.1


shock_tube["p"] = shock_tube["p"] * ureg("barye")
shock_tube["rho"] = shock_tube["rho"] * ureg("g/cc")
write_initial_hdf5(filename="initial_conditions", initial_condition_dict=shock_tube)

# Plot the results
try:
    fig, (ax1, ax2) = plt.subplots(figsize=(18, 8), nrows=2, ncols=1)
    for ax, v in zip([ax1, ax2], ["rho", "p"]):
        vc = ax.plot(shock_tube["xc"][:, 1], shock_tube[v][:, 1], "-o")
        ax.set_ylabel(v)
        ax.set_xlabel("X")
    plt.show()

    # plt.figure(figsize=(10, 2))
    # plt.pcolormesh(
    #     shock_tube["x"].m,
    #     shock_tube["y"].m,
    #     shock_tube["rho"].m,
    #     ec="k",
    #     lw=0.1,
    #     antialiased=True,
    #     cmap="viridis",
    # )
    # plt.axis("equal")
    # plt.show()
except Exception:
    pass
