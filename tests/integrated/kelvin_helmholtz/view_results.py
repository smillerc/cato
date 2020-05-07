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
from pycato import load_2d_dataset

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
ds = load_2d_dataset(".")

plt.figure(figsize=(12, 12))
ds.density[-1].plot(x="x", y="y")
t = ds.time[-1].data
plt.title(
    f"Kelvin-Helmholtz Test @ {now} \nsimulation t={t:.2f} s \nwalltime={walltime_sec} s\nbranch: {branch} \ncommit: {short_hash}"
)

plt.axis("equal")
plt.tight_layout()
plt.savefig("kh_2d_results.png")
