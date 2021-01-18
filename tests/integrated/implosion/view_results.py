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

fig, (contour_ax, resid_ax) = plt.subplots(ncols=2, nrows=1, figsize=(24, 12))
ds.density[-1].plot.pcolormesh(
    x="x",
    y="y",
    ec="k",
    lw=0.1,
    antialiased=True,
    cmap="viridis",
    # vmin=0.0,
    # vmax=2.4e-3,
    ax=contour_ax,
)

ds.density[-1].plot.contour(
    x="x", y="y", colors="k", linewidths=0.5, antialiased=True, levels=12, ax=contour_ax
)

# Plot the residual history
df.plot(
    x="time",
    y=["rho", "rho_u", "rho_v", "rho_E"],
    kind="line",
    logy=True,
    ax=resid_ax,
    antialiased=True,
    lw=0.75,
)
resid_ax.set_ylabel("residual")
resid_ax.set_xlabel("time")
resid_ax.set_ylim(1e-16, 0.1)

t = ds.time[-1].data
contour_ax.set_title(
    f"Implosion Test @ {now} \nsimulation t={t:.4f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash} \nscheme: {scheme}"
)

contour_ax.axis("equal")
plt.tight_layout()
plt.savefig("implosion_2d_results.png")