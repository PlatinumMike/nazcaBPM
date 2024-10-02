#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 19:38:06 2024

Loading and plotting the field data

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt


def get_field(filename: str, Nx: int, Ny: int, cmpx: bool = True):
    data = np.genfromtxt(filename, delimiter=",")
    if cmpx:
        data = data[:, 0] + 1.0j * data[:, 1]
    data = np.reshape(data, shape=(Nx, Ny))
    return data


nx = 40
ny = 60
nz = 4728
Lx = 4
Ly = 6
Lz = 40
x = np.linspace(-0.5 * Lx, 0.5 * Lx, nx)
y = np.linspace(-0.5 * Ly, 0.5 * Ly, ny)
z = np.linspace(0, Lz, nz)

initial_field = get_field("../build/" + "initial_field.csv", Nx=nx, Ny=ny)
final_field = get_field("../build/" + "final_field.csv", Nx=nx, Ny=ny)
field_slice = get_field("../build/" + "field_slice.csv", Nx=nz, Ny=ny)

index_start = get_field("../build/" + "index_start.csv", Nx=nx, Ny=ny, cmpx=False)
index_end = get_field("../build/" + "index_end.csv", Nx=nx, Ny=ny, cmpx=False)
index_cross = get_field("../build/" + "index_cross.csv", Nx=ny, Ny=nz, cmpx=False)


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
plt.contourf(z, y, intensity.T, 20, cmap="inferno")
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
