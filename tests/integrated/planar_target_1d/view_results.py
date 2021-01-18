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
from pint import UnitRegistry

ureg = UnitRegistry()

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

# try:
scheme = f"{ds.attrs['scheme_flux_solver']}({ds.attrs['scheme_spatial_reconstruction']} {ds.attrs['scheme_limiter']})"
# except Exception:
#     scheme = None

fig, ax = plt.subplots(figsize=(12, 6))


ds["mach_x"] = ds.x_velocity / ds.sound_speed
ds.density[-1].plot(x="x", label="Density")
# ds.pressure[-1].plot(x="x", label="Pressure")

plt.legend()
ax2 = ax.twinx()
ds.pressure[-1].plot(x="x", label="Pressure", ax=ax2, color="k")

ax3 = ax.twinx()
ax3.spines["right"].set_position(("axes", 1.1))
# ds.mach_x[-1].plot(x="x", label="Mach", ax=ax3, color="r")
ds.x_velocity[-1].plot(x="x", label="Velocity", ax=ax3, color="r")

ax3.set_title(None)
ax2.set_title(None)
ax.set_title(None)

time = ds.time[-1].time.data * ureg(ds.time_units)
plt.title(f"time={time.to('ps'):.2f~}")
plt.tight_layout()
plt.savefig("layered_target_1d_results.png")
