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
from datetime import datetime

sys.path.append("../../..")
from pycato import *

pp_dir = "./post_processed_results"
if not os.path.exists(pp_dir):
    os.mkdir(pp_dir)
# Load cato results
ds = load_dataset(".", use_dask=True)

try:
    scheme = f"{ds.attrs['scheme.flux_solver']} + {ds.attrs['scheme.limiter']}"
except Exception:
    scheme = None

t_range = np.arange(0, 4.1, 0.1)

for t in t_range:
    fig, (dens_ax, mach_ax) = plt.subplots(ncols=1, nrows=2, figsize=(24, 12))
    ds.density.sel(time=t, method="nearest").plot.pcolormesh(
        x="x", y="y", antialiased=True, cmap="viridis", ax=dens_ax,
    )

    ds.density.sel(time=t, method="nearest").plot.contour(
        x="x",
        y="y",
        colors="k",
        linewidths=0.5,
        antialiased=True,
        levels=12,
        ax=dens_ax,
    )

    ds.mach_y.sel(time=t, method="nearest").plot.pcolormesh(
        x="x", y="y", antialiased=True, cmap="RdBu_r", ax=mach_ax,
    )

    ds.mach_y.sel(time=t, method="nearest").plot.contour(
        x="x",
        y="y",
        colors="k",
        linewidths=0.5,
        antialiased=True,
        levels=12,
        ax=mach_ax,
    )

    now = datetime.now()
    current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    dens_ax.set_title(f"Path: {os.path.abspath(os.curdir)}, Time: {current_time}")
    plt.tight_layout()
    plt.savefig(f"{pp_dir}/results_t={t:.3f}ns.png")
    plt.close()
