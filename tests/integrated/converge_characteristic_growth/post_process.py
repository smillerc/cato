# -*- coding: utf-8 -*-
import h5py
import numpy as np
import os
import pytz
from datetime import datetime
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd

sys.path.append("../../..")
from pycato import *


perturbed_ds = load_dataset(folder="perturbed")
symmetric_ds = load_dataset(folder="symmetric")

diff = perturbed_ds.x_velocity[-1] - symmetric_ds.x_velocity[-1]

diff.plot()
plt.show()
