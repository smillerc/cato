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

tz = pytz.timezone("America/New_York")
now = datetime.now(tz)

try:
    short_hash = (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("utf-8")
        .strip()
    )
except Exception:
    short_hash = "N/A"

try:
    branch = (
        subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"])
        .decode("utf-8")
        .strip()
    )

    if branch == "HEAD":
        branch = os.getenv("CI_COMMIT_REF_NAME")
    if branch is None:
        branch = "N/A"

except Exception:
    branch = "N/A"

try:
    walltime_sec = np.loadtxt("timing.csv", delimiter=",", skiprows=1)[-1][1]
except Exception:
    walltime_sec = "N/A"

# Load cato results
from dask.diagnostics import ProgressBar

with ProgressBar():
    ds = load_multiple_steps("results/step*.h5", ini_file="input.ini")

try:
    scheme = f"{ds.attrs['scheme_flux_solver']}({ds.attrs['scheme_spatial_reconstruction']} {ds.attrs['scheme_limiter']})"
except Exception:
    scheme = None

df = pd.read_csv("residual_hist.csv", index_col=False)

ds["mach_y"] = ds.sound_speed / ds.y_velocity
ds.mach_y[-1].plot.pcolormesh(
    x="x", y="y", ec="k", lw=0.1, antialiased=True, vmin=-1e-3, vmax=1e-3, cmap='RdBu_r'
)

ds.mach_y[-1].plot.contour(x="x", y="y", colors="k", linewidths=0.5, antialiased=True)

t = ds.time[-1].data
plt.title(
    f"Sedov Test @ {now} \nsimulation t={t:.4f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash} \nscheme: {scheme}"
)

plt.axis("equal")

try:
    plt.show()
except Exception:
    pass
