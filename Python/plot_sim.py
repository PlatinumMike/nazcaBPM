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


def get_index(filename: str, grid1: str, grid2: str):
    f = h5py.File(filename, "r")
    index = f["refractive_index"][:]
    _grid1 = f[grid1][:]
    _grid2 = f[grid2][:]
    data = {"index": index, grid1: _grid1, grid2: _grid2}
    f.close()
    return data


def get_field(filename: str, grid1: str, grid2: str):
    f = h5py.File(filename, "r")
    real = f["real"][:]
    imag = f["imag"][:]
    _grid1 = f[grid1][:]
    _grid2 = f[grid2][:]
    f.close()
    data = {"field": real + 1.0j * imag, grid1: _grid1, grid2: _grid2}
    return data


# %% load files
field_start = get_field("../build/field_start.h5", "xgrid", "ygrid")
field_end = get_field("../build/field_end.h5", "xgrid", "ygrid")
field_yz = get_field("../build/field_yz.h5", "ygrid", "zgrid")
field_xz = get_field("../build/field_xz.h5", "xgrid", "zgrid")

index_start = get_index("../build/index_start.h5", "xgrid", "ygrid")
index_end = get_index("../build/index_end.h5", "xgrid", "ygrid")
index_xz = get_index("../build/index_xz.h5", "xgrid", "zgrid")
index_yz = get_index("../build/index_yz.h5", "ygrid", "zgrid")

# %% plotting

plt.figure()
plt.contourf(field_start["ygrid"], field_start["xgrid"], field_start["field"].real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

plt.figure()
plt.contourf(field_end["ygrid"], field_end["xgrid"], field_end["field"].real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

intensity_yz = np.abs(field_yz["field"]) ** 2
intensity_xz = np.abs(field_xz["field"]) ** 2

# warning, the y direction is mirrored, because ymin is at the top, ymax at the bottom.
plt.figure()
plt.contourf(field_yz["zgrid"], field_yz["ygrid"], intensity_yz, 20, cmap="inferno")
plt.xlabel("z")
plt.ylabel("y")
plt.show()

plt.figure()
plt.contourf(field_xz["zgrid"], field_xz["xgrid"], intensity_xz, 20, cmap="inferno")
plt.xlabel("z")
plt.ylabel("x")
plt.show()


plt.figure()
plt.contourf(index_end["ygrid"], index_end["xgrid"], index_end["index"])
plt.xlabel("y")
plt.ylabel("x")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_start["ygrid"], index_start["xgrid"], index_start["index"])
plt.xlabel("y")
plt.ylabel("x")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_yz["zgrid"], index_yz["ygrid"], index_yz["index"])
plt.xlabel("z")
plt.ylabel("y")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_xz["zgrid"], index_xz["xgrid"], index_xz["index"])
plt.xlabel("z")
plt.ylabel("x")
plt.colorbar()
plt.show()
