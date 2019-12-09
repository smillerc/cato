# -*- coding: utf-8 -*-
"""Methods for generating initial grids for fvleg_2d
in the .h5 format
"""
import numpy as np
import h5py


def make_uniform_grid(n_nodes, xrange, yrange):
    """Generate a uniform grid. This will output a dictionary
    that contains the appropriate arrays, which include the ghost
    cell layer.

    Parameters
    ----------
    n_nodes : tuple
        Number of non-ghost nodes in each direction (x,y)
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

    x = np.linspace(xrange[0], xrange[1], n_nodes[0])
    ldx = x[1] - x[0]
    rdx = x[-1] - x[-2]
    x = np.array([x[0] - ldx] + list(x) + [x[-1] + rdx])

    y = np.linspace(yrange[0], yrange[1], n_nodes[1])
    ldy = y[1] - y[0]
    rdy = y[-1] - y[-2]
    y = np.array([y[0] - ldy] + list(y) + [y[-1] + rdy])

    xc = np.zeros(x.shape[0] - 1)
    yc = np.zeros(y.shape[0] - 1)

    xx, yy = np.meshgrid(x, y)
    xxc, yyc = np.meshgrid(xc, yc)

    rho = np.ones_like(xxc)
    pressure = np.ones_like(xxc)
    u = np.ones_like(xxc)
    v = np.ones_like(xxc)

    dy = (np.diff(yy[:, 0]) / 2)[0]
    dx = (np.diff(xx[0, :]) / 2)[0]
    xc = xx[:-1, :-1] + dx
    yc = yy[:-1, :-1] + dy

    return {
        "x": xx,
        "y": yy,
        "rho": rho,
        "u": u,
        "v": v,
        "p": pressure,
        "xc": xc,
        "yc": yc,
    }


def write_initial_hdf5(filename, initial_condition_dict, boundary_conditions_dict):
    """ Write the initial conditions to an hdf5 file.

    Parameters
    ----------
    filename : str
        Name of the grid file
    initial_condition_dict : dictionary
        Dictionary created by the `make_uniform_grid` method. Must
        contain the following keys ['x', 'y', 'rho', 'u', 'v', 'p']
    boundary_conditions_dict : dictionary
        Dictionary containing the boundary conditions
        at ['+x','+y','-x','-y']
    """

    if not filename.endswith(".h5"):
        filename += ".h5"

    print("Writing to: ", filename)
    with h5py.File(filename, mode="w") as h5:
        h5.create_dataset(
            "/x",
            data=initial_condition_dict["x"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/y",
            data=initial_condition_dict["y"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/density",
            data=initial_condition_dict["rho"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/x_velocity",
            data=initial_condition_dict["u"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/y_velocity",
            data=initial_condition_dict["v"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/pressure",
            data=initial_condition_dict["p"].astype(np.float32),
            compression="gzip",
            compression_opts=9,
        )

        # Boundary Conditions
        h5.create_dataset("/plus_x_bc", data=boundary_conditions_dict["+x"])
        h5.create_dataset("/plus_y_bc", data=boundary_conditions_dict["+y"])
        h5.create_dataset("/minus_x_bc", data=boundary_conditions_dict["-x"])
        h5.create_dataset("/minus_y_bc", data=boundary_conditions_dict["-y"])
        h5.close()
