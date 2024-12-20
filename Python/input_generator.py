#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 20:48:56 2024

Generate JSON input files for the BPM code.
You can create your own json files from scratch, but it is recommended to use this class because then the input
gets validated by pydantic automatically. 
This way you are less likely to waste CPU time simulating some incorrect settings.
(fail-fast principle)

@author: mike
"""

from typing import List, Tuple

from pydantic import BaseModel


class Layer(BaseModel):
    name: str
    zmin: float
    zmax: float
    index: float


class XS(BaseModel):
    xs_name: str = "default"
    layer_list: List[Layer] = []
    background_index: float = 1.0

    def get_index(self, z: float) -> float:
        for layer in self.layer_list:
            if z >= layer.zmin and z <= layer.zmax:
                return layer.index
        return self.background_index


class Port(BaseModel):
    name: str = "a0"
    placement: str = "left"
    yspan: float = 4.0
    zspan: float = 4.0
    y0: float = 0.0
    z0: float = 0.0
    port_resolution_y: int = 30
    port_resolution_z: int = 30


class Shape(BaseModel):
    cell_name: str = "cell"
    poly: List[Tuple[float, float]] = []
    xs_name: str = "default"


class SettingsBPM(BaseModel):
    reference_index: float
    wl: float
    resolution_x: int
    resolution_y: int
    resolution_z: int
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float
    shapes: List[Shape] = []
    input_ports: List[Port]
    output_ports: List[Port]
    cross_sections: List[XS]
    absolute_path_output: str = "/mnt/workdir"
    scheme_parameter: float = 0.5
    pml_strength: float = 5.0
    pml_thickness: float = 1.0
    dry_run: bool = False
