# -*- coding: utf-8 -*-
"""
A simple script to view the results from the simulation
"""

import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import os, sys

sys.path.append("../../..")
from pycato import load_1d_dataset

t = 0.17

# Load cato results
ds = load_1d_dataset(".")

plt.figure(figsize=(12, 6))

ds.density.sel(time=t, method="nearest").plot()
plt.title(f"Time: {t:.2f}")
plt.ylabel("Density")
plt.xlabel("X")
plt.savefig("shu_osher_1d_results.png")
