# -*- coding: utf-8 -*-
import numpy as np
import h5py


def make_uniform_grid(n_nodes=(1, 1), xrange=(0, 1), yrange=(0, 1)):
    """Create an empty grid"""
    ghost_cell = 1

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
    """
    boundary_conditions_dict: Dict
        Needs to have the keys [+x,+y,-x,-y]
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


bc_dict = {"+x": "periodic", "+y": "periodic", "-x": "periodic", "-y": "periodic"}

# Make a simple grid with dummy initial conditions
simple = make_uniform_grid(n_nodes=(5, 3), xrange=(-2, 2), yrange=(-2, 2))
simple["rho"] = simple["rho"] * 1
simple["u"] = simple["p"] * -2
simple["v"] = simple["p"] * 3
simple["p"] = simple["p"] * 4

write_initial_hdf5(
    filename="simple", initial_condition_dict=simple, boundary_conditions_dict=bc_dict
)
