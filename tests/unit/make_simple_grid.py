# -*- coding: utf-8 -*-
import numpy as np
import h5py
import pint

ureg = pint.UnitRegistry()


def make_uniform_grid(n_cells, xrange, yrange, n_ghost_layers=2):
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
simple = make_uniform_grid(n_cells=(4, 4), xrange=(-4, 4), yrange=(-2, 2), n_ghost_layers=2)


# !   Domain cartoon layout
# Note: this is what the unit tests will expect (w/o the ghost layer population)
# ! j  |---|---||---|---|---|---||---|---|
# ! 6  | 0 | 6 || 8 | 0 | 0 | 6 || 8 | 0 | <- j=6; aka top_ghost (grid%ubounds_halo(2))
# !    |---|---||---|---|---|---||---|---|
# ! 5  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 |
# !    |===|===||===|===|===|===||===|===|
# ! 4  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 | <- j=4; aka top (grid%ubounds(2))
# !    |---|---||---|---|---|---||---|---|
# ! 3  | 0 | 6 || 8 | 0 | 0 | 6 || 8 | 0 |
# !    |---|---||---|---|---|---||---|---|
# ! 2  | 0 | 6 || 8 | 0 | 0 | 6 || 8 | 0 |
# !    |---|---||---|---|---|---||---|---|
# ! 1  | 5 | 2 || 1 | 5 | 5 | 2 || 1 | 5 | <- j=1; aka bottom (grid%lbounds(2))
# !    |===|===||===|===|===|===||===|===|
# ! 0  | 7 | 3 || 4 | 7 | 7 | 3 || 4 | 7 |
# !    |---|---||---|---|---|---||---|---|
# ! -1 | 0 | 6 || 8 | 0 | 0 | 6 || 8 | 0 | <- j=-1; aka bottom_ghost (grid%lbounds_halo(2))
# !    |---|---||---|---|---|---||---|---|
# !     -1   0    1   2   3   4    5   6    i

data = np.zeros_like(simple["rho"])

data[2:6,2] = [1, 5, 5, 2]
data[2:6,3] = [8, 9, 9, 6]
data[2:6,4] = [8, 9, 9, 6]
data[2:6,5] = [4, 7, 7, 3]

print(data.T)
simple["rho"] = data * ureg('g/cc')
simple["u"] = -data * ureg('cm/s')
simple["v"] = data * ureg('cm/s')
simple["p"] = data * ureg('barye')

write_initial_hdf5(filename="simple", initial_condition_dict=simple)
