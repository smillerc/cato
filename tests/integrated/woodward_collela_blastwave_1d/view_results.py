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

# Remove the ghost layers (bc's)
ds = ds.where(ds["ghost_cell"] == 0, drop=True)
serialize_dataset(ds)

t = 0.17
actual_time = ds.density.sel(t=t, method="nearest").t.data

plt.figure(figsize=(12, 6))
ds.density.sel(t=t, method="nearest").plot(x="x")

plt.title(
    f"Shu-Osher 1D Test @ {now} \nsimulation t={actual_time:.2f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash}"
)
plt.ylabel("Density [g/cc]")
plt.xlabel("X [cm]")
plt.legend()
plt.ylim(0, 6)
plt.tight_layout()
plt.savefig("shu_osher_1d_results.png")

try:
    plt.show()
except Exception:
    pass
