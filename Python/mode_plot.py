#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 21:21:04 2024

Load and plot field of mode solver

@author: mike
"""

import os

import numpy as np
import matplotlib.pyplot as plt
import h5py


def get_field(filename: str, grid1: str, grid2: str):
    f = h5py.File(filename, "r")
    real = f["real"][:]
    imag = f["imag"][:]
    _grid1 = f[grid1][:]
    _grid2 = f[grid2][:]
    f.close()
    data = {"field": real + 1.0j * imag, grid1: _grid1, grid2: _grid2}
    return data


def plot_slice(field: dict):
    plt.figure()
    plt.contourf(field["ygrid"], field["zgrid"], np.abs(field["field"].T) ** 2, 20)
    plt.xlabel("y")
    plt.ylabel("z")
    plt.title("|u|^2")
    plt.colorbar()
    plt.show()


dirname = os.path.dirname(__file__)
build_dir = os.path.join(os.path.dirname(dirname), "simulations")

port = "a0"
iterations = 61
fields = []
for i in range(iterations):
    field_xy = get_field(os.path.join(build_dir, f"monitor_{port}_yz_{i}.h5"), "ygrid", "zgrid")
    fields.append(field_xy)


plot_slice(fields[-1])
