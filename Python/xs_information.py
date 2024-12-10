#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 21:11:59 2024

Cross section class

@author: mike
"""
from typing import List


class Layer:
    def __init__(
        self,
        name: str,
        zmin: float,
        zmax: float,
        index: float,
    ):
        self.name = name
        self.zmin = zmin
        self.zmax = zmax
        self.index = index

    def layer2str(self):
        return {"name": self.name, "zmin": self.zmin, "zmax": self.zmax, "index": self.index}


class XS:
    def __init__(self, name: str, layer_list: List[Layer], background_index: float = 1.0):
        self.xs_name = name
        self.layer_list = layer_list
        self.background_index = background_index

    def xs2str(self):
        layer_list_str = [layer.layer2str() for layer in self.layer_list]
        return {
            "xs_name": self.xs_name,
            "layer_list": layer_list_str,
            "background_index": self.background_index,
        }


# some examples below. The "default" XS just a cross section with no core layers. This will be used in the simulation for points that lie outside of the polygons.
def soi_strip(height: float = 0.220):
    core = Layer("core", -height / 2, height / 2, 3.5)
    xs_core = XS("soi-strip", [core], 1.5)
    xs_default = XS("default", [], 1.5)
    return xs_core, xs_default


def sin_strip(height: float = 0.350):
    core = Layer("core", -height / 2, height / 2, 2.0)
    xs_core = XS("sin-strip", [core], 1.5)
    xs_default = XS("default", [], 1.5)
    return xs_core, xs_default


def InP_rib(height: float = 0.400, height_base: float = 1.0, height_top: float = 1.5):
    index_SiO2 = 1.50
    index_core = 3.40
    index_clad = 3.16
    height_sub = 1.0e6  # just something very large so it extends through the entire bbox.

    substrate = Layer("substrate", -height_sub, 0.0, index_clad)
    z0 = 0.0
    nInp = Layer("nInp", z0, height_base, index_clad)
    z0 += height_base
    core = Layer("core", z0, z0 + height, index_core)
    z0 += height
    pInp = Layer("pInp", z0, z0 + height_top, index_clad)
    xs_core = XS("sin-strip", [substrate, nInp, core, pInp], index_SiO2)
    xs_default = XS("default", [substrate], index_SiO2)
    return xs_core, xs_default
