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
import seaborn as sns

cp = sns.color_palette()

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

ds = ds.squeeze().drop("y").where(ds.x < 10.0, drop=True)
ds['u+cs'] = ds.x_velocity + ds.sound_speed
# try:
scheme = f"{ds.attrs['scheme_flux_solver']}({ds.attrs['scheme_spatial_reconstruction']} {ds.attrs['scheme_limiter']})"
# except Exception:
#     scheme = None

t = 9.25
actual_time = ds.density.sel(time=t, method="nearest").time.data


fig, ax = plt.subplots(figsize=(10,10))
ds.density.sel(time=t, method='nearest').plot(x='x', ax=ax, label='Density', color=cp[0])

ax.set_ylabel('Density')
ax.set_xlabel('X')
ax.set_title(None)

# Pressure
ax3 = ax.twinx()
ds.pressure.sel(time=t, method='nearest').plot(x='x', ax=ax3, color=cp[2])
ax3.spines["right"].set_position(("axes", 1.1))
ax3.set_ylabel('Pressure', color=cp[2])
ax3.set_title(None)

# Velocity
ax2 = ax.twinx()
ds['u+cs'].sel(time=t, method='nearest').plot(x='x', ax=ax2, color=cp[1])
ax2.set_ylabel("U+Cs", color=cp[1])

plt.title(
    f"Sod 1D Test @ {now} \nsimulation t={actual_time:.2f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash} \nscheme: {scheme}"
)

plt.tight_layout()
plt.savefig("kink_test.png")

try:
    plt.show()
except Exception:
    pass
