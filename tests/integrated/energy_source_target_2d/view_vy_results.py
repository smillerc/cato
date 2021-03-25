# -*- coding: utf-8 -*-
"""
A simple script to view the results from the simulation
"""

import h5py
import numpy as np
import os
import pytz
from datetime import datetime
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd

sys.path.append("../../..")
from pycato import *

# Load cato results
from dask.diagnostics import ProgressBar

with ProgressBar():
    ds = load_multiple_steps("results/step*.h5", ini_file="input.ini")


vmin = -0.1
vmax = 0.1
lvls = 64
plt.figure(figsize=(16,4))
ds.y_velocity.sel(time=0.46, method="nearest").plot.contourf(
    levels=lvls,
    x="x",
    y="y",
    antialiased=True,
    vmin=vmin,
    vmax=vmax,
    cmap="RdBu_r",
)

ds.y_velocity.sel(time=0.46, method="nearest").plot.contour(
    levels=lvls,
    x="x",
    y="y",
    vmin=vmin,
    vmax=vmax,
    colors="k",
    linewidths=0.5,
    antialiased=True,
)
plt.xlim(40, 62)
plt.tight_layout()
plt.savefig('vy.png')
plt.show()
