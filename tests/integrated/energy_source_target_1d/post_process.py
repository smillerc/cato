# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
from dask.diagnostics import ProgressBar

sys.path.append("../../..")
from pycato import *

pp_dir = Path.cwd() / "post_processed_results"

if not pp_dir.exists():
    pp_dir.mkdir()

with ProgressBar():
    perturbed = load_multiple_steps(f"perturbed/step*.h5", use_dask=True,)
with ProgressBar():
    symmetric = load_multiple_steps(f"symmetric/step*.h5", use_dask=True,)

xmin = 65
xmax = 79
tmin = 0.1

pert_v = (
    perturbed.x_velocity.squeeze()
    .drop("y")
    .where(
        (perturbed.x > xmin) & (perturbed.x < xmax) & (perturbed.time > tmin),
        drop=True,
    )
)
sym_v = (
    symmetric.x_velocity.squeeze()
    .drop("y")
    .where(
        (symmetric.x > xmin) & (symmetric.x < xmax) & (symmetric.time > tmin),
        drop=True,
    )
)

dv = sym_v - pert_v

fig, ax = plt.subplots(figsize=(16, 8))
ax2 = ax.twinx()
for t in np.arange(0.3, 4.1, 0.1):
    pert_v.sel(time=t, method="nearest").plot(x="x", label=f"Time: {t:.3f}", ax=ax)
    dv.sel(time=t, method="nearest").plot(x="x", label=f"Time: {t:.3f}", ax=ax2)
plt.legend()
plt.savefig(f"{pp_dir}/slices.png")
plt.close()

plt.figure(figsize=(16, 8))
plt.title("Max delta v vs time")
dv.max(dim="i").plot()
plt.savefig(f"{pp_dir}/final.png")
