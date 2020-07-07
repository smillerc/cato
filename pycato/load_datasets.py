# -*- coding: utf-8 -*-
"""Read the .h5 files from a simulation into an xarray Dataset"""
import configparser
import os

import h5py
import numpy as np
import pint
import xarray as xr


from .unit_registry import ureg


def as_dict(config):
    """
    Converts a ConfigParser object into a dictionary.

    The resulting dictionary has sections as keys which point to a dict of the
    sections options as key => value pairs.
    """
    the_dict = {}
    for section in config.sections():
        the_dict[section] = {}
        for key, val in config.items(section):
            try:
                the_dict[section][key] = float(val.strip("'"))
            except ValueError:
                the_dict[section][key] = val.strip("'")

    return the_dict


def flatten_dict(dd, separator=".", prefix=""):
    """Flatten a dictionary"""
    return (
        {
            prefix + separator + k if prefix else k: v
            for kk, vv in dd.items()
            for k, v in flatten_dict(vv, separator, kk).items()
        }
        if isinstance(dd, dict)
        else {prefix: dd}
    )


def read_ini(filename):
    """Read the input .ini file

    Parameters
    ----------
    filename : str
        Path of the input.ini file

    Returns
    -------
    Dict
        A flattened dictionary of the input options
    """
    from configparser import ConfigParser

    config = ConfigParser()
    config.read(filename)
    return flatten_dict(as_dict(config))


def load_dataset(folder, use_dask=False, drop_ghost=True):
    """Read the dataset info

    Parameters
    ----------
    folder : [type]
        [description]
    """

    results_folder = os.path.join(folder, "results")

    # Read the input file
    input_dict = read_ini(os.path.join(folder, "input.ini"))

    step_files = sorted(
        [
            os.path.join(results_folder, f)
            for f in os.listdir(results_folder)
            if f.endswith(".h5") and f.startswith("step")
        ]
    )
    h5_files = [h5py.File(f, "r") for f in step_files]
    h5_files

    var_list = [
        # "ghost_cell",
        "density",
        "pressure",
        "sound_speed",
        # "temperature",
        "x_velocity",
        "y_velocity",
    ]

    data_vars = {}
    space_dims = ("t", "x", "y")

    # The ghost cell array is a special case
    data_vars[f"ghost_cell"] = xr.Variable(
        ("x", "y"),
        np.array(h5_files[0]["/ghost_cell"][()].T, dtype=np.int8),
        attrs={"units": h5_files[0]["/ghost_cell"].attrs["units"].decode("utf-8")},
    )

    if use_dask:
        import dask.array as da

        for v in var_list:
            chunk_size = h5_files[0][f"/{v}"].shape
            datasets = [h5[f"/{v}"][()] for h5 in h5_files]
            arrays = [da.from_array(dataset, chunks=chunk_size) for dataset in datasets]
            data_array = da.transpose(da.stack(arrays, axis=0), axes=[0, 2, 1])

            data_vars[f"{v}"] = xr.Variable(
                space_dims,
                data_array,
                attrs={"units": h5_files[0][f"/{v}"].attrs["units"].decode("utf-8")},
            )

    else:

        for v in var_list:
            data_vars[f"{v}"] = xr.Variable(
                space_dims,
                np.array([h5[f"/{v}"][()].T for h5 in h5_files], dtype=np.float32),
                attrs={"units": h5_files[0][f"/{v}"].attrs["units"].decode("utf-8")},
            )

    x = np.array(h5_files[0][f"/x"][()].T, dtype=np.float32)
    x_units = h5_files[0][f"/x"].attrs["units"].decode("utf-8")
    y = np.array(h5_files[0][f"/y"][()].T, dtype=np.float32)

    # Get the cell centers
    dy = (np.diff(x[0, :]) / 2.0)[0]
    dx = (np.diff(y[:, 0]) / 2.0)[0]

    # cell center locations
    xc = x[:-1, 0] + dx
    yc = y[0, :-1] + dy

    coords = {
        "t": np.array([h5[f"/time"][()] for h5 in h5_files], dtype=np.float32),
        "x": xc,
        "y": yc,
    }
    time_units = h5_files[0][f"/time"].attrs["units"].decode("utf-8")

    # Get the details about the CATO build
    info_attr = {}
    info = [
        "build_type",
        "compile_hostname",
        "compile_os",
        "compiler_flags",
        "compiler_version",
        "git_changes",
        "git_hash",
        "git_ref",
        "version",
    ]
    for v in info:
        try:
            info_attr[v] = h5_files[0]["/cato_info"].attrs[f"{v}"].decode("utf-8")
        except Exception:
            pass

    attr_dict = info_attr
    attr_dict.update(input_dict)
    attr_dict["time_units"] = time_units
    attr_dict["x_units"] = x_units

    ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attr_dict)

    # Compute a few variables
    ds["mach_x"] = ds.x_velocity / ds.sound_speed
    ds["mach_y"] = ds.y_velocity / ds.sound_speed

    if drop_ghost:
        return ds.where(ds["ghost_cell"] == 0, drop=True)
    else:
        return ds


def serialize_dataset(dataset):
    """Serialize the xarray dataset to file

    Parameters
    ----------
    dataset : xr.Dataset
    """
    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in dataset.data_vars}
    filename = "dataset.nc"
    print(f"Saving dataset to: {os.path.abspath(filename)}")
    dataset.to_netcdf(filename, engine="h5netcdf", encoding=encoding)
