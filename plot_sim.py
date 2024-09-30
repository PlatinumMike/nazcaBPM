#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 19:38:06 2024

Loading and plotting the field data

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt


def get_field(filename, Nx, Ny):
    data = np.genfromtxt(filename, delimiter=",")
    data = data[:, 0] + 1.0j * data[:, 1]
    data = np.reshape(data, shape=(Nx, Ny))
    return data


initial_field = get_field("build/" + "initial_field.csv", Nx=120, Ny=180)
final_field = get_field("build/" + "final_field.csv", Nx=120, Ny=180)
field_slice = get_field("build/" + "field_slice.csv", Nx=109839, Ny=180)

plt.figure()
plt.contourf(initial_field.real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

plt.figure()
plt.contourf(final_field.real, 20)
plt.xlabel("y")
plt.ylabel("x")
plt.show()

intensity = np.abs(field_slice) ** 2

# warning, the y direction is mirrored, because ymin is at the top, ymax at the bottom.
plt.figure()
plt.contourf(intensity.T, 20, cmap="inferno")
plt.xlabel("z")
plt.ylabel("y")
plt.show()
