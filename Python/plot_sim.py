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
field_yz_start = get_field("../build/field_yz_start.h5", "ygrid", "zgrid")
field_yz_end = get_field("../build/field_yz_end.h5", "ygrid", "zgrid")
field_xz = get_field("../build/field_xz.h5", "xgrid", "zgrid")
field_xy = get_field("../build/field_xy.h5", "xgrid", "ygrid")

index_yz_start = get_index("../build/index_yz_start.h5", "ygrid", "zgrid")
index_yz_end = get_index("../build/index_yz_end.h5", "ygrid", "zgrid")
index_xy = get_index("../build/index_xy.h5", "xgrid", "ygrid")
index_xz = get_index("../build/index_xz.h5", "xgrid", "zgrid")

xend = field_xz["xgrid"][-1]

# %% plotting

plt.figure()
plt.contourf(field_yz_start["ygrid"], field_yz_start["zgrid"], field_yz_start["field"].real.T, 20)
plt.xlabel("y")
plt.ylabel("z")
plt.title("Re(u), x=0")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(field_yz_end["ygrid"], field_yz_end["zgrid"], field_yz_end["field"].real.T, 20)
plt.xlabel("y")
plt.ylabel("z")
plt.title(f"Re(u), x={xend}")
plt.colorbar()
plt.show()

intensity_xz = np.abs(field_xz["field"]) ** 2
intensity_xy = np.abs(field_xy["field"]) ** 2

plt.figure()
plt.contourf(field_xz["xgrid"], field_xz["zgrid"], intensity_xz.T, 20, cmap="inferno")
plt.xlabel("x")
plt.ylabel("z")
plt.title("Intensity")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(field_xy["xgrid"], field_xy["ygrid"], intensity_xy.T, 20, cmap="inferno")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Intensity")
plt.colorbar()
plt.show()


plt.figure()
plt.contourf(index_yz_end["ygrid"], index_yz_end["zgrid"], index_yz_end["index"].T)
plt.xlabel("y")
plt.ylabel("z")
plt.title("Refractive index")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_yz_start["ygrid"], index_yz_start["zgrid"], index_yz_start["index"].T)
plt.xlabel("y")
plt.ylabel("z")
plt.title("Refractive index")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_xz["xgrid"], index_xz["zgrid"], index_xz["index"].T)
plt.xlabel("x")
plt.ylabel("z")
plt.title("Refractive index")
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(index_xy["xgrid"], index_xy["ygrid"], index_xy["index"].T)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Refractive index")
plt.colorbar()
plt.show()
