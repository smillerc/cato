# -*- coding: utf-8 -*-
"""Read the .h5 files from a simulation into an xarray Dataset"""
import h5py
import numpy as np
import os
import xarray as xr
import pint

ureg = pint.UnitRegistry()


def read_stepfile(file):
    """Read a single .hdf5 file

    Parameters
    ----------
    file : string
        Filename of the hdf5 file

    Returns
    -------
    Dict
        Dictionary of the data from within the hdf5 file
    """
    data = {}

    var_list = ["x", "density", "x_velocity", "y_velocity", "pressure"]

    with h5py.File(file, "r") as h5:
        # Transpose to match the index convention within the code
        data["time"] = h5["/time"][()]
        for var in var_list:
            try:
                data[var] = h5[f"/{var}"][()].T[:, 1]
            except Exception:
                print(f"Unable to read {var}")
    return data


def read_1d_dataset(folder, units="cgs"):
    """[summary]

    Parameters
    ----------
    folder : [type]
        [description]
    units : str, optional
        [description], by default 'cgs'

    Returns
    -------
    [type]
        [description]

    Raises
    ------
    Exception
        [description]
    """

    if units != "cgs" and units != "icf":
        raise Exception("Unsupported units not in the set = {cgs, icf}")

    step_files = sorted(
        [
            os.path.join(folder, f)
            for f in os.listdir(folder)
            if f.endswith(".h5") and f.startswith("step")
        ]
    )
    data = read_stepfile(step_files[0])
    data_dim_x = data["density"].shape[0]
    data_dim_t = len(step_files)

    data = {
        "x": np.zeros((data_dim_t, data_dim_x + 1)),
        "time": np.zeros(data_dim_t),
        "density": np.zeros((data_dim_t, data_dim_x)),
        "pressure": np.zeros((data_dim_t, data_dim_x)),
        "x_velocity": np.zeros((data_dim_t, data_dim_x)),
        "y_velocity": np.zeros((data_dim_t, data_dim_x)),
    }

    for t, f in enumerate(step_files):
        single_step_data = read_stepfile(f)
        data["x"][t, :] = single_step_data["x"]
        data["time"][t] = single_step_data["time"]
        data["x_velocity"][t, :] = single_step_data["x_velocity"]
        data["y_velocity"][t, :] = single_step_data["y_velocity"]
        data["density"][t, :] = single_step_data["density"]
        data["pressure"][t, :] = single_step_data["pressure"]

    density = data["density"] * ureg("g/cc")
    pressure = data["pressure"] * ureg("barye")
    x_vel = data["x_velocity"] * ureg("cm/s")
    y_vel = data["y_velocity"] * ureg("cm/s")
    time = data["time"] * ureg("s")
    xc = np.cumsum(np.diff(data["x"][0])) * ureg("cm")

    if units == "cgs":
        density_units = "g/cc"
        pressure_units = "barye"
        vel_units = "cm/s"
        time_units = "s"
        length_units = "cm"

    elif units == "icf":
        density_units = "g/cc"
        pressure_units = "Mbar"
        vel_units = "um/ns"
        time_units = "ns"
        length_units = "um"

    ds = xr.Dataset(
        data_vars={
            "Density": xr.Variable(
                ("time", "x"),
                density.to(density_units).m,
                attrs={"units": density_units},
            ),
            "Pressure": xr.Variable(
                ("time", "x"),
                pressure.to(pressure_units).m,
                attrs={"units": pressure_units},
            ),
            "x_velocity": xr.Variable(
                ("time", "x"), x_vel.to(vel_units).m, attrs={"units": vel_units}
            ),
            "y_velocity": xr.Variable(
                ("time", "x"), y_vel.to(vel_units).m, attrs={"units": vel_units}
            ),
        },
        coords={"time": time.to(time_units).m, "x": xc.to(length_units).m},
    )
    return ds
