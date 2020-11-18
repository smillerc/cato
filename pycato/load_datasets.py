# -*- coding: utf-8 -*-
"""Read the .h5 files from a simulation into an xarray Dataset"""
import configparser
import os

import h5py
import numpy as np
import pint
import xarray as xr
import dask
import dask.array as da
from glob import glob
from pathlib import Path


from .unit_registry import ureg

var_dict = {
    "density": {"long_name": "Density", "standard_name": "rho"},
    "pressure": {"long_name": "Pressure", "standard_name": "p"},
    "x_velocity": {"long_name": "X-Velocity", "standard_name": "u"},
    "y_velocity": {"long_name": "Y-Velocity", "standard_name": "v"},
    "sound_speed": {"long_name": "Sound Speed", "standard_name": "cs"},
}


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


def flatten_dict(dd, separator="_", prefix=""):
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
    flattened_dict = flatten_dict(as_dict(config))
    for k, v in flattened_dict.items():
        if isinstance(v, str):
            flattened_dict[k] = v.replace("'", "").replace('"', "")
    return flattened_dict


class _MultiFileCloser:
    __slots__ = ("file_objs",)

    def __init__(self, file_objs):
        self.file_objs = file_objs

    def close(self):
        for f in self.file_objs:
            f.close()


def load_single(file, drop_ghost=True, use_dask=True, var_list="all", ini_file=None):
    """Load a single step file and generate an xarray Dataset

    Parameters
    ----------
    file : str or Path
        Location of the file to load
    drop_ghost : bool, optional
        Drop all of the ghost cells, by default True
    var_list : List, optional
        Load only a specific set of variables, by default 'all'

    Returns
    -------
    xarray Dataset
    """

    if var_list == "all":
        var_list = [
            "density",
            "pressure",
            "sound_speed",
            "x_velocity",
            "y_velocity",
            "ghost_cell",
        ]

    data_vars = {}
    space_dims = ("i", "j")

    if not file.endswith(".h5"):
        raise Exception("Step files must be .h5 files")

    h5 = h5py.File(file, "r")

    for v in var_list:
        try:
            h5[f"/{v}"].shape
        except KeyError:
            continue

        if use_dask:
            chunk_size = h5[f"/{v}"].shape
            array = da.from_array(h5[f"/{v}"], chunks=chunk_size)
            array = da.transpose(array)
        else:
            array = h5[f"/{v}"][()].T.astype(np.float32)

        try:
            long_name = var_dict[v]["long_name"]
        except Exception:
            long_name = ""

        try:
            description = h5[f"/{v}"].attrs["description"].decode("utf-8")
        except Exception:
            description = ""

        try:
            standard_name = var_dict[v]["standard_name"]
        except Exception:
            standard_name = ""

        try:
            units = h5[f"/{v}"].attrs["units"].decode("utf-8")
        except Exception:
            units = ""

        data_vars[f"{v}"] = xr.Variable(
            space_dims,
            array,
            attrs={
                "units": units,
                "description": description,
                "long_name": long_name,
                "standard_name": standard_name,
            },
        )

    x = h5[f"/x"][()].T.astype(np.float32)
    x_units = h5[f"/x"].attrs["units"].decode("utf-8")
    y = h5[f"/y"][()].T.astype(np.float32)

    # Get the cell centers
    dy = (np.diff(x[0, :]) / 2.0)[0]
    dx = (np.diff(y[:, 0]) / 2.0)[0]

    # cell center locations
    xc = x[:-1, 0] + dx
    yc = y[0, :-1] + dy

    coords = {
        "time": h5[f"/time"][()].astype(np.float32),
        "x": (["i"], xc),
        "y": (["j"], yc),
    }

    time_units = h5[f"/time"].attrs["units"].decode("utf-8")

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
            info_attr[v] = h5["/cato_info"].attrs[f"{v}"].decode("utf-8")
        except Exception:
            pass

    attr_dict = info_attr
    attr_dict["time_units"] = time_units
    attr_dict["space_units"] = x_units

    if ini_file:
        input_dict = read_ini(ini_file)
        attr_dict.update(input_dict)

    ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attr_dict)
    if ini_file:
        try:
            ds.attrs["title"] = ds.attrs["general_title"]
        except Exception:
            pass

    if drop_ghost:
        try:
            ds = ds.where(ds["ghost_cell"] == 0, drop=True)
            return ds.drop("ghost_cell")
        except KeyError:
            return ds
    else:
        return ds


def load_multiple_steps(paths, use_dask=True, **kwargs):
    """Load multiple datasets and concatenate them together. This is very similar
    to xarray's open_mfdataset

    Parameters
    ----------
    paths : str or List
        Either a string glob in the form ``"path/to/my/files/*.h5"`` or an explicit list of
        files to open. Paths can be given as strings or as pathlib Paths.

    Returns
    -------
    xarray.Dataset
    """

    if isinstance(paths, str):
        paths = sorted(glob(paths))
    else:
        paths = [str(p) if isinstance(p, Path) else p for p in paths]

    if use_dask:
        # Use dask's delayed function on loading single step files
        open_ = dask.delayed(load_single)
        getattr_ = dask.delayed(getattr)
        datasets = [open_(f, **kwargs) for f in paths]

        file_objs = [getattr_(ds, "_file_obj") for ds in datasets]
        datasets, file_objs = dask.compute(datasets, file_objs)
    else:
        datasets = [load_single(f, use_dask=False, **kwargs) for f in paths]

    # Concatenate all of the datasets together based on the time dimension
    combined = xr.concat(datasets, dim="time")

    if use_dask:
        combined._file_obj = _MultiFileCloser(file_objs)
    return combined


def serialize_compressed_dataset(dataset, filename="dataset.nc", lvl=6):
    """Serialize the xarray dataset to file

    Parameters
    ----------
    dataset : xr.Dataset
    filename : str
    lvl : int
        Level of compression to use
    """
    comp = dict(zlib=True, complevel=lvl)
    encoding = {var: comp for var in dataset.data_vars}
    print(f"Saving dataset to: {os.path.abspath(filename)}")
    dataset.to_netcdf(filename, engine="h5netcdf", encoding=encoding)
