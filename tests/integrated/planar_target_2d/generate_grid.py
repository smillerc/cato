import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pint
import h5py

ureg = pint.UnitRegistry()


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


def linear_spacing(start, initial_dx, max_dx, scale, end, reverse=False):
    dx = initial_dx
    x = [start]
    while True:
        new_x = x[-1] + dx
        x.append(new_x)
        dx = dx * scale
        dx = min(dx, max_dx)
        if x[-1] > end:
            break

    x_new = np.array(x, dtype=np.float64)
    x_scale = x_new / x_new[-1] * end

    x_out = x_scale + start
    if reverse:
        x = [start] + list(np.cumsum(-np.diff(np.flip(x))))
    dx_last = x[-1] - x[-2]
    dx_first = x[1] - x[0]
    ncells = len(x)
    return x, dx_first, dx_last, ncells


def x_spacing(left_size, fine_size, right_size, min_dx, max_dx, scale=1.09):

    x_left, dx_first, dx_last, ncells = linear_spacing(
        start=0,
        initial_dx=min_dx,
        scale=scale,
        max_dx=max_dx,
        end=left_size,
        reverse=True,
    )
    x_middle = np.arange(x_left[-1], x_left[-1] + fine_size, dx_last)
    x_right, dx_first, dx_last, ncells = linear_spacing(
        start=x_middle[-1],
        initial_dx=dx_last,
        max_dx=max_dx,
        scale=scale,
        end=x_middle[-1] + right_size,
    )

    x_mesh = np.unique(np.concatenate([x_left, x_middle, x_right]))

    return x_mesh


def add_halo(mesh, n_halo):
    ldx = mesh[1] - mesh[0]
    rdx = mesh[-1] - mesh[-2]

    if n_halo == 1:
        new_mesh = np.array(
            [mesh[0] - ldx] + list(mesh) + [mesh[-1] + rdx], dtype=np.float64
        )
    elif n_halo == 2:
        new_mesh = np.array(
            [mesh[0] - 2 * ldx]
            + [mesh[0] - ldx]
            + list(mesh)
            + [mesh[-1] + rdx]
            + [mesh[-1] + 2 * rdx],
            dtype=np.float64,
        )
    elif n_halo == 3:
        new_mesh = np.array(
            [mesh[0] - 3 * ldx]
            + [mesh[0] - 2 * ldx]
            + [mesh[0] - ldx]
            + list(mesh)
            + [mesh[-1] + rdx]
            + [mesh[-1] + 2 * rdx]
            + [mesh[-1] + 3 * rdx],
            dtype=np.float64,
        )
    else:
        raise Exception("Unable to work with n_halo that isn't 1 2, or 3")

    return new_mesh


# Main inputs
# ---------------------------------------------
plot = True
gamma = 5.0 / 3.0
n_halo = 2


layers_thickness = [70, 35, 8, 55]

initial_pressure = 1e9  # barye

gas_density = 0.001  # g/cc
gas_pressure = initial_pressure

ice_density = 0.25  # g/cc
ice_pressure = initial_pressure

ch_density = 1.0  # g/cc
ch_pressure = initial_pressure

vac_density = 0.01  # g/cc
vac_pressure = 40e12  # barye
apply_vac_pressure_gradient = True
vac_pressure_end = 1e10
apply_vac_vel = False  # initialize the corona w/ an outflow velocity


# Ice layer perturbation
apply_pert = True
x_pert = [85]
lambda_x = 2
lambda_y = 20
ymax = lambda_y / 2.0

min_dx = 0.1
max_dx = 1e4
dy = ymax / 1
# dy = 0.1

perturbation_frac = 0.5

# ---------------------------------------------

spacing = [30, 60 + 10 + 8,  60]

if np.sum(spacing) != np.sum(layers_thickness):
    raise Exception(
        "Double check spacing and layer thicknesses, they don't equal the same amount"
    )

x_mesh = x_spacing(
    left_size=spacing[0],
    fine_size=spacing[1],
    right_size=spacing[2],
    min_dx=min_dx,
    max_dx=max_dx,
    scale=1.5,
)

x_mesh = add_halo(x_mesh, n_halo)
y_mesh = np.linspace(0, ymax, num=int(ymax / dy) + 1, endpoint=True)
y_mesh = add_halo(y_mesh, n_halo)

xc = x_mesh[:-1] + np.diff(x_mesh) / 2.0
yc = y_mesh[:-1] + np.diff(y_mesh) / 2.0

nj = yc.shape[0]
ni = xc.shape[0]

rho = np.ones((ni, nj))
u = np.ones((ni, nj))
v = np.ones((ni, nj))
p = np.ones((ni, nj))
pert = np.ones((ni, nj))

coords = {"x": (["i"], xc), "y": (["j"], yc)}
ds = xr.Dataset(
    data_vars=dict(
        density=(["i", "j"], rho),
        x_velocity=(["i", "j"], u),
        y_velocity=(["i", "j"], v),
        pressure=(["i", "j"], p),
        perturbation=(["i", "j"], pert),
    ),
    coords=coords,
)

layer_x = np.cumsum(layers_thickness)
gas_layer = ds.x < layer_x[0]
ice_layer = (ds.x >= layer_x[0]) & (ds.x <= layer_x[1])
ch_layer = (ds.x >= layer_x[1]) & (ds.x <= layer_x[2])
vac_layer = ds.x > layer_x[2]

ds.density[gas_layer] = gas_density
ds.density[ice_layer] = ice_density
ds.density[ch_layer] = ch_density
ds.density[vac_layer] = vac_density

ds.pressure.data[:, :] = gas_pressure
ds.pressure.data[vac_layer] = vac_pressure

# add a gradient in the coronal pressure
if apply_vac_pressure_gradient:
    p_s = vac_pressure
    p_e = vac_pressure_end
    
    cx = ds.x[vac_layer]

    cx_len = (cx - cx[0])
    cx_last = cx[-(n_halo+1)]
    cx_max = cx[-1]
    cx[:] = cx_len / cx_max
    slope = (p_e - p_s) / cx_last
    p_c = slope * cx + p_s

    ds.pressure.data[vac_layer] = np.tile(p_c, nj).reshape(nj, len(p_c)).T
    for j in range(n_halo+1):
        ds.pressure.data[-j,:] = vac_pressure_end

ds.x_velocity.data[:, :] = 0
ds.y_velocity.data[:, :] = 0

if apply_vac_vel:
    v_shell = np.sqrt(2.0 / (gamma + 1.0) * vac_pressure / ch_density)
    ds.x_velocity[vac_layer] = -v_shell


k_x = (2 * np.pi) / lambda_x
k_y = (2 * np.pi) / lambda_y

Xn, Yn = np.meshgrid(x_mesh, y_mesh)
Xn = Xn.T
Yn = Yn.T

Xc, Yc = np.meshgrid(xc, yc)
Xc = Xc.T
Yc = Yc.T

if apply_pert:
    for x0 in x_pert:
        pert = perturbation_frac * np.abs(
            np.exp(-k_x * ((Xc - x0)) ** 2) * -((1.0 - np.cos(k_y * Yc)) / 2.0)
        )

        pert[pert < 1e-6] = 0.0
        ds.perturbation.data[:, :] = pert
        ds.perturbation[np.logical_not(ice_layer)] = 0.0
        print(f"Maximum perturbation fraction: {pert.max():.3f}")

        # Apply to the density
        ds.density[ice_layer] = ds.density[ice_layer] * (1 - ds.perturbation[ice_layer])

print(f"Total # Cells: {ds.density.size:,}")
if plot:
    # fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(16, 10))
    # ds.density.plot.pcolormesh(
    #     x="x", y="y", lw=0.8, alpha=0.05, antialiased=True, cmap="RdBu_r", ax=ax1,
    # )
    # ds.perturbation.plot.pcolormesh(
    #     x="x", y="y", lw=0.8, alpha=0.1, antialiased=True, cmap="RdBu_r", ax=ax2,
    # )
    # ax1.axis("equal")
    # ax2.axis("equal")
    # ax1.set_ylim(0, ymax)
    # ax2.set_ylim(0, ymax)
    # plt.show()

    fig, ax = plt.subplots(figsize=(16, 10))
    plt.plot(x_mesh[:-1], np.diff(x_mesh), color="k")
    plt.xlabel("X [um]")
    plt.ylabel(r"$\delta x$ (mesh spacing)")
    plt.yscale("log")
    ax2 = ax.twinx()
    ds.density[:, -3].plot(x="x", ax=ax2, marker='s', ms=3)
    # ds.pressure[:, -3].plot(x="x", ax=ax2, marker='s', ms=3)
    ax2.set_ylabel("Density [g/cc]")
    ax2.set_title(None)
    ax.set_title("Mesh Spacing")
    plt.show()

initial_condition_dict = {
    "x": (Xn * ureg("um")).to("cm"),
    "y": (Yn * ureg("um")).to("cm"),
    "rho": ds.density.data * ureg("g/cc"),
    "u": ds.x_velocity.data * ureg("cm/s"),
    "v": ds.y_velocity.data * ureg("cm/s"),
    "p": ds.pressure.data * ureg("barye"),
    "n_ghost_layers": n_halo,
}

write_initial_hdf5("initial_conditions", initial_condition_dict)
