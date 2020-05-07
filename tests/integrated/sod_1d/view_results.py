# -*- coding: utf-8 -*-
"""
A simple script to view the results from the simulation
"""

import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import sod

sys.path.append("../../..")
from pycato import load_1d_dataset

t = 0.2
gamma = 1.4
npts = 500

# exact results
positions, regions, values = sod.solve(
    left_state=(1, 1, 0),
    right_state=(0.1, 0.125, 0.0),
    geometry=(0.0, 1.0, 0.5),
    t=t,
    gamma=gamma,
    npts=npts,
)
p = values["p"]
rho = values["rho"]
u = values["u"]

# Load cato results
ds = load_1d_dataset(".")

plt.figure(figsize=(12, 6))

ds.density.sel(time=t, method="nearest").plot(label="CATO Density")
plt.plot(values["x"], rho, label="Exact Density")

ds.velocity.sel(time=t, method="nearest").plot(label="CATO Velocity")
plt.plot(values["x"], u, label="Exact Veclocity")

ds.pressure.sel(time=t, method="nearest").plot(label="CATO Pressure")
plt.plot(values["x"], p, label="Exact Pressure")

plt.title(f"Time: {t:.2f}")
plt.ylabel("")
plt.xlabel("X")

plt.legend()
plt.savefig("sod_1d_results.png")
