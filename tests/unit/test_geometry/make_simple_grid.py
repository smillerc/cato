# -*- coding: utf-8 -*-
import numpy as np
import h5py
import pint

ureg = pint.UnitRegistry()


def make_uniform_grid(n_cells, xrange, yrange, n_ghost_layers=1):
    """Generate a uniform grid. This will output a dictionary
    that contains the appropriate arrays, which include the ghost
    cell layer.

    Parameters
    ----------
    n_cells : tuple
        Number of non-ghost cells in each direction (x,y)
    xrange : tuple
        Extent of the domain in x (min,max)
    yrange : tuple
        Extent of the domain in y (min,max)

    Returns
    -------
    dictionary
        A dictionary that contains the conserved variables (rho, u velocity, v
        velocity, p), grid (x,y) points, and the cell center (xc,yc) points
    """

    dx = float(xrange[1] - xrange[0]) / float(n_cells[0])
    dy = float(yrange[1] - yrange[0]) / float(n_cells[1])

    x = np.linspace(
        start=xrange[0] - n_ghost_layers * dx,
        stop=xrange[1] + n_ghost_layers * dx,
        num=(n_cells[0] + 1) + (n_ghost_layers * 2),
        dtype=np.float64,
    )

    y = np.linspace(
        start=yrange[0] - n_ghost_layers * dy,
        stop=yrange[1] + n_ghost_layers * dy,
        num=(n_cells[1] + 1) + (n_ghost_layers * 2),
        dtype=np.float64,
    )

    xc = np.zeros(x.shape[0] - 1, dtype=np.float64)
    yc = np.zeros(y.shape[0] - 1, dtype=np.float64)

    # 2d versions
    y_2d, x_2d = np.meshgrid(y, x)  # nodes

    # cell-centered arrays
    # node_shape = (x_2d.shape[0], x_2d.shape[1])
    cell_shape = (x_2d.shape[0] - 1, x_2d.shape[1] - 1)
    rho = np.ones(cell_shape)
    u = np.ones(cell_shape)
    v = np.ones(cell_shape)
    p = np.ones(cell_shape)

    # cell center locations
    xc = x_2d[:-1, :-1] + dx
    yc = y_2d[:-1, :-1] + dy

    return {
        "x": x_2d * ureg("cm"),
        "y": y_2d * ureg("cm"),
        "rho": rho * ureg("g/cc"),
        "u": u * ureg("cm/s"),
        "v": v * ureg("cm/s"),
        "p": p * ureg("barye"),
        "xc": xc * ureg("cm"),
        "yc": yc * ureg("cm"),
        "n_ghost_layers": n_ghost_layers,
    }


def write_initial_hdf5(filename, initial_condition_dict):
    """ Write the initial conditions to an hdf5 file.

    Parameters
    ----------
    filename : str
        Name of the grid file
    initial_condition_dict : dictionary
        Dictionary created by the `make_uniform_grid` method. Must
        contain the following keys ['x', 'y', 'rho', 'u', 'v', 'p']
    """

    if not filename.endswith(".h5") or not filename.endswith(".hdf5"):
        filename += ".h5"

    print("Writing to: ", filename)
    with h5py.File(filename, mode="w") as h5:

        data = initial_condition_dict["n_ghost_layers"]
        h5.create_dataset("/n_ghost_layers", data=data)

        data = initial_condition_dict["x"].to("cm").m
        dset = h5.create_dataset(
            "/x", data=data.astype(np.float64).T, compression="gzip", compression_opts=9
        )
        dset.attrs["units"] = "cm"

        data = initial_condition_dict["y"].to("cm").m
        dset = h5.create_dataset(
            "/y", data=data.astype(np.float64).T, compression="gzip", compression_opts=9
        )
        dset.attrs["units"] = "cm"

        data = initial_condition_dict["rho"].to("g/cc").m
        dset = h5.create_dataset(
            "/density",
            data=data.astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        dset.attrs["units"] = "g/cc"

        data = initial_condition_dict["u"].to("cm/s").m
        dset = h5.create_dataset(
            "/x_velocity",
            data=data.astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        dset.attrs["units"] = "cm/s"

        data = initial_condition_dict["v"].to("cm/s").m
        dset = h5.create_dataset(
            "/y_velocity",
            data=data.astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        dset.attrs["units"] = "cm/s"

        data = initial_condition_dict["p"].to("barye").m
        dset = h5.create_dataset(
            "/pressure",
            data=data.astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        dset.attrs["units"] = "barye"

        h5.close()


# Make a simple grid with dummy initial conditions
simple = make_uniform_grid(n_cells=(4, 2), xrange=(-2, 2), yrange=(-2, 2))
simple["rho"] = simple["rho"] * 1
simple["u"] = simple["u"] * -2
simple["v"] = simple["v"] * 3
simple["p"] = simple["p"] * 4

write_initial_hdf5(filename="simple", initial_condition_dict=simple)
