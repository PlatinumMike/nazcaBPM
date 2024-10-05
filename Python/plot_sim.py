#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 19:38:06 2024

Loading and plotting the field data

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py


def get_index(filename: str):
    f = h5py.File(filename, "r")
    index = f["refractive_index"][:]
    f.close()
    return index


def get_field(filename: str):
    f = h5py.File(filename, "r")
    real = f["real"][:]
    imag = f["imag"][:]
    f.close()
    return real + 1.0j * imag


# %% load files

filename = "../build/index_data.h5"


f = h5py.File(filename, "r")

x = f["xgrid"][:]
y = f["ygrid"][:]
z = f["zgrid"][:]
index3d = f["refractive_index"][:]

f.close()


index_start = get_index("../build/index_start.h5")
index_end = get_index("../build/index_end.h5")
index_cross = get_index("../build/index_cross.h5")
field_slice = get_field("../build/field_slice.h5")
final_field = get_field("../build/final_field.h5")
initial_field = get_field("../build/initial_field.h5")

# %% plotting

plt.figure()
plt.contourf(y, x, initial_field.real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

plt.figure()
plt.contourf(y, x, final_field.real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

intensity = np.abs(field_slice) ** 2

# warning, the y direction is mirrored, because ymin is at the top, ymax at the bottom.
plt.figure()
plt.contourf(z, y, intensity, 20, cmap="inferno")
plt.xlabel("z")
plt.ylabel("y")
plt.show()


plt.figure()
plt.contourf(y, x, index_end)
plt.xlabel("y")
plt.ylabel("x")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(z, y, index_cross)
plt.xlabel("z")
plt.ylabel("y")
plt.colorbar()
plt.show()
