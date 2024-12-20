#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 21:11:59 2024

Cross section class

@author: mike
"""
from typing import List
from input_generator import Layer, XS


# some examples below. The "default" XS just a cross section with no core layers. This will be used in the simulation for points that lie outside of the polygons.
def soi_strip(height: float = 0.220):
    core = Layer(name="core", zmin=-height / 2, zmax=height / 2, index=3.5)
    xs_core = XS(xs_name="soi-strip", layer_list=[core], background_index=1.5)
    xs_default = XS(xs_name="default", layer_list=[], background_index=1.5)
    return xs_core, xs_default


def sin_strip(height: float = 0.350):
    core = Layer(name="core", zmin=-height / 2, zmax=height / 2, index=2.0)
    xs_core = XS(xs_name="sin-strip", layer_list=[core], background_index=1.5)
    xs_default = XS(xs_name="default", layer_list=[], background_index=1.5)
    return xs_core, xs_default


def InP_rib(height: float = 0.400, height_base: float = 1.0, height_top: float = 1.5):
    index_SiO2 = 1.50
    index_core = 3.40
    index_clad = 3.16
    height_sub = 1.0e6  # just something very large so it extends through the entire bbox.

    substrate = Layer(name="substrate", zmin=-height_sub, zmax=0.0, index=index_clad)
    z0 = 0.0
    nInp = Layer(name="nInp", zmin=z0, zmax=height_base, index=index_clad)
    z0 += height_base
    core = Layer(name="core", zmin=z0, zmax=z0 + height, index=index_core)
    z0 += height
    pInp = Layer(name="pInp", zmin=z0, zmax=z0 + height_top, index=index_clad)
    xs_core = XS(xs_name="sin-strip", layer_list=[substrate, nInp, core, pInp], background_index=index_SiO2)
    xs_default = XS(xs_name="default", layer_list=[substrate], background_index=index_SiO2)
    return xs_core, xs_default


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    zgrid = np.linspace(-2.0, 5.0, 100)
    xs_core, xs_default = InP_rib(height=0.35)
    index_profile = [xs_core.get_index(z) for z in zgrid]

    plt.figure()
    plt.plot(zgrid, index_profile)
    plt.xlabel("z (um)")
    plt.ylabel("Refractive index")
    plt.show()
