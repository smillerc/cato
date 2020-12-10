import h5py
import numpy as np


data = np.ones((8,8))
with h5py.File("data.h5", mode="w") as h5:
    dset = h5.create_dataset("/density",data=data)
    dset.attrs["units"] = "g/cc"
    dset.attrs["description"] = "Mass Density"