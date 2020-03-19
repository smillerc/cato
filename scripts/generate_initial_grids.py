# -*- coding: utf-8 -*-
"""Methods for generating initial grids for cato
in the .h5 format
"""
import numpy as np
import h5py
import pint

ureg = pint.UnitRegistry()


def make_uniform_grid(n_cells, xrange, yrange):
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

    x = np.linspace(xrange[0], xrange[1], n_cells[0] - 1, dtype=np.float64)
    ldx = x[1] - x[0]
    rdx = x[-1] - x[-2]
    x = np.array([x[0] - ldx] + list(x) + [x[-1] + rdx], dtype=np.float64)

    y = np.linspace(yrange[0], yrange[1], n_cells[1] - 1, dtype=np.float64)
    ldy = y[1] - y[0]
    rdy = y[-1] - y[-2]
    y = np.array([y[0] - ldy] + list(y) + [y[-1] + rdy], dtype=np.float64)

    xc = np.zeros(x.shape[0] - 1, dtype=np.float64)
    yc = np.zeros(y.shape[0] - 1, dtype=np.float64)

    xx, yy = np.meshgrid(x, y)
    xxc, yyc = np.meshgrid(xc, yc)

    rho = np.ones_like(xxc, dtype=np.float64)
    pressure = np.ones_like(xxc, dtype=np.float64)
    u = np.ones_like(xxc, dtype=np.float64)
    v = np.ones_like(xxc, dtype=np.float64)

    dy = (np.diff(yy[:, 0]) / 2.0)[0]
    dx = (np.diff(xx[0, :]) / 2.0)[0]
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


def make_2d_layered_grid(
    layer_thicknesses,
    layer_n_cells,
    layer_density,
    layer_u,
    layer_v,
    layer_pressure,
    y_thickness,
):
    """Create a 2d layered grid (uniform in j), layers are in i

    Parameters
    ----------
    layer_thicknesses : 1d array
        Array with units for the thickness of each layer
    layer_n_cells : list
        List with the number of cells for each layer
    layer_density : 1d array
        Array with units for the density of each layer
    layer_u : 1d array
        Array with units for the x velocity of each layer
    layer_v : 1d array
        Array with units for the y velocity of each layer
    layer_pressure : 1d array
        Array with units for the pressure of each layer

    Returns
    -------
    dictionary
        A dictionary that contains the conserved variables (rho, u velocity, v
        velocity, p), grid points (x, y), and the cell center (xc, yc) points
    """

    cumulative_thickness = 0
    cumulative_ncells = 1  # start at 1 due to the ghost cell
    node_x = []
    layer_cell_idx_ranges = []  # [(first cell, last cell), etc]

    layer_thicknesses = layer_thicknesses.to("cm").m
    for ncells, thickness in zip(layer_n_cells, layer_thicknesses):
        start_x = cumulative_thickness
        stop_x = start_x + thickness
        cumulative_thickness += thickness

        start = cumulative_ncells
        stop = start + ncells
        layer_cell_idx_ranges.append((start, stop - 1))
        cumulative_ncells += ncells

        layer_x = np.linspace(start_x, stop_x, ncells + 1, endpoint=True)
        node_x.append(layer_x)

    # X location of all of the nodes/vertices
    x = np.concatenate([xl[:-1] for xl in node_x] + [[node_x[-1][-1]]])

    # right/left cell size - for making ghost layers
    ldx = x[1] - x[0]
    rdx = x[-1] - x[-2]

    # find the minimum cell spacing
    dx = np.diff(x).min()
    n_ycells = int(y_thickness / dx)

    y = np.linspace(0, y_thickness, n_ycells, dtype=np.float64)
    ldy = y[1] - y[0]
    rdy = y[-1] - y[-2]
    y = np.array([y[0] - ldy] + list(y) + [y[-1] + rdy], dtype=np.float64)

    # new x location for nodes/vertices including ghost cells
    x = np.array([x[0] - ldx] + list(x) + [x[-1] + rdx], dtype=np.float64)

    # 2d versions
    y_2d, x_2d = np.meshgrid(y, x)  # nodes

    # cell-centered arrays
    # node_shape = (x_2d.shape[0], x_2d.shape[1])
    cell_shape = (x_2d.shape[0] - 1, x_2d.shape[1] - 1)
    rho = np.zeros(cell_shape)
    u = np.zeros(cell_shape)
    v = np.zeros(cell_shape)
    p = np.zeros(cell_shape)

    # Assign rho, u, v, and p to the cell-centered arrays
    for l, layer in enumerate(layer_cell_idx_ranges):
        s = layer[0]
        e = layer[1] + 1

        # Convert to cgs units
        rho[s:e, :] = layer_density[l].to("g/cc").m
        u[s:e, :] = layer_u[l].to("cm/s").m
        v[s:e, :] = layer_v[l].to("cm/s").m
        p[s:e, :] = layer_pressure[l].to("barye").m

    # Assign ghost layer values
    for i in [0, -1]:
        rho[i, :] = layer_density[i].to("g/cc").m
        u[i, :] = layer_u[i].to("cm/s").m
        v[i, :] = layer_v[i].to("cm/s").m
        p[i, :] = layer_pressure[i].to("barye").m

    # cell spacing
    dy = (np.diff(y_2d[0, :]) / 2.0)[0]
    dx = (np.diff(x_2d[:, 0]) / 2.0)[0]

    # cell center locations
    xc = x_2d[:-1, :-1] + dx
    yc = y_2d[:-1, :-1] + dy

    return {
        "x": x_2d,
        "y": y_2d,
        "rho": rho,
        "u": u,
        "v": v,
        "p": p,
        "xc": xc,
        "yc": yc,
    }


def make_1d_layered_grid(
    layer_thicknesses, layer_n_cells, layer_density, layer_u, layer_v, layer_pressure
):
    """Create a 1d layered grid

    Parameters
    ----------
    layer_thicknesses : 1d array
        Array with units for the thickness of each layer
    layer_n_cells : list
        List with the number of cells for each layer
    layer_density : 1d array
        Array with units for the density of each layer
    layer_u : 1d array
        Array with units for the x velocity of each layer
    layer_v : 1d array
        Array with units for the y velocity of each layer
    layer_pressure : 1d array
        Array with units for the pressure of each layer

    Returns
    -------
    dictionary
        A dictionary that contains the conserved variables (rho, u velocity, v
        velocity, p), grid points (x, y), and the cell center (xc, yc) points
    """

    cumulative_thickness = 0
    cumulative_ncells = 1  # start at 1 due to the ghost cell
    node_x = []
    layer_cell_idx_ranges = []  # [(first cell, last cell), etc]

    layer_thicknesses = layer_thicknesses.to("cm").m
    for ncells, thickness in zip(layer_n_cells, layer_thicknesses):
        start_x = cumulative_thickness
        stop_x = start_x + thickness
        cumulative_thickness += thickness

        start = cumulative_ncells
        stop = start + ncells
        layer_cell_idx_ranges.append((start, stop - 1))
        cumulative_ncells += ncells

        layer_x = np.linspace(start_x, stop_x, ncells + 1, endpoint=True)
        node_x.append(layer_x)

    # X location of all of the nodes/vertices
    x = np.concatenate([xl[:-1] for xl in node_x] + [[node_x[-1][-1]]])

    # right/left cell size - for making ghost layers
    ldx = x[1] - x[0]
    rdx = x[-1] - x[-2]

    # Since this is a "1d" grid, it needs 3 cells in y
    # (1 on either side for ghost layers)
    y = np.array([-ldx, 0, ldx, ldx * 2], dtype=np.float64) - ldx / 2

    # new x location for nodes/vertices including ghost cells
    x = np.array([x[0] - ldx] + list(x) + [x[-1] + rdx], dtype=np.float64)

    # 2d versions
    y_2d, x_2d = np.meshgrid(y, x)  # nodes

    # cell-centered arrays
    rho = np.zeros((x.shape[0] - 1, 3))
    u = np.zeros((x.shape[0] - 1, 3))
    v = np.zeros((x.shape[0] - 1, 3))
    p = np.zeros((x.shape[0] - 1, 3))

    # Assign rho, u, v, and p to the cell-centered arrays
    for l, layer in enumerate(layer_cell_idx_ranges):
        s = layer[0]
        e = layer[1] + 1

        # Convert to cgs units
        rho[s:e, :] = layer_density[l].to("g/cc").m
        u[s:e, :] = layer_u[l].to("cm/s").m
        v[s:e, :] = layer_v[l].to("cm/s").m
        p[s:e, :] = layer_pressure[l].to("barye").m

    # Assign ghost layer values
    for i in [0, -1]:
        rho[i, :] = layer_density[i].to("g/cc").m
        u[i, :] = layer_u[i].to("cm/s").m
        v[i, :] = layer_v[i].to("cm/s").m
        p[i, :] = layer_pressure[i].to("barye").m

    # cell spacing
    dy = (np.diff(y_2d[0, :]) / 2.0)[0]
    dx = (np.diff(x_2d[:, 0]) / 2.0)[0]

    # cell center locations
    xc = x_2d[:-1, :-1] + dx
    yc = y_2d[:-1, :-1] + dy

    return {
        "x": x_2d,
        "y": y_2d,
        "rho": rho,
        "u": u,
        "v": v,
        "p": p,
        "xc": xc,
        "yc": yc,
    }


def make_1d_in_x_uniform_grid(n_cells, limits=(0, 1)):
    """Generate a uniform grid 1d grid in x. This will output a dictionary
    that contains the appropriate arrays, which include the ghost
    cell layer.

    Parameters
    ----------
    n_cells : tuple
        Number of non-ghost cells
    limits : tuple
        Extent of the domain (min,max)

    Returns
    -------
    dictionary
        A dictionary that contains the conserved variables (rho, u velocity, v
        velocity, p), grid points, and the cell center (xc, yc) points
    """

    x = np.linspace(limits[0], limits[1], n_cells - 1, dtype=np.float64)
    ldx = x[1] - x[0]
    rdx = x[-1] - x[-2]
    x = np.array([x[0] - ldx] + list(x) + [x[-1] + rdx], dtype=np.float64)

    y = np.array([-ldx, 0, ldx, ldx * 2], dtype=np.float64) - ldx / 2

    # make empty arrays for the cell center locations
    xc = np.zeros(x.shape[0] - 1, dtype=np.float64)
    yc = np.zeros(y.shape[0] - 1, dtype=np.float64)

    # 2d versions
    xx, yy = np.meshgrid(x, y)  # nodes
    xxc, yyc = np.meshgrid(xc, yc)  # cell centers

    rho = np.ones_like(xxc, dtype=np.float64)
    pressure = np.ones_like(xxc, dtype=np.float64)
    u = np.ones_like(xxc, dtype=np.float64)
    v = np.ones_like(xxc, dtype=np.float64)

    # cell spacing
    dy = (np.diff(yy[:, 0]) / 2.0)[0]
    dx = (np.diff(xx[0, :]) / 2.0)[0]

    # cell center locations
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

    if not filename.endswith(".h5"):
        filename += ".h5"

    print("Writing to: ", filename)
    with h5py.File(filename, mode="w") as h5:
        h5.create_dataset(
            "/x",
            data=initial_condition_dict["x"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/y",
            data=initial_condition_dict["y"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/density",
            data=initial_condition_dict["rho"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/x_velocity",
            data=initial_condition_dict["u"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/y_velocity",
            data=initial_condition_dict["v"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "/pressure",
            data=initial_condition_dict["p"].astype(np.float64).T,
            compression="gzip",
            compression_opts=9,
        )

        h5.close()
