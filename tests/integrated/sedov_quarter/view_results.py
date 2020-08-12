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

# # Load cato results
from dask.diagnostics import ProgressBar

with ProgressBar():
    ds = load_multiple_steps("results/step*.h5", ini_file="input.ini")


plt.figure(figsize=(12, 12))
ds.density[-1].plot.pcolormesh(
    x="x",
    y="y",
    ec="k",
    lw=0.1,
    antialiased=True,
    cmap="viridis",
    vmin=0.0,
    vmax=2.4e-3,
)
t = ds.time[-1].data
plt.title(
    f"Sedov Test @ {now} \nsimulation t={t:.4f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash}"
)

# Plot a circle to check symmetry with the eye
theta = np.linspace(0, np.pi / 2.0, 50)
radius = 0.131313
circle_x = radius * np.cos(theta)
circle_y = radius * np.sin(theta)
plt.plot(circle_x, circle_y, color="k")

plt.axis("equal")
r = 0.2
plt.xlim(0, r)
plt.ylim(0, r)
plt.tight_layout()
plt.savefig("sedov_2d_results.png")

try:
    plt.show()
except Exception:
    pass
