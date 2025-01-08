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

from enum import Enum
from typing import List, Tuple, Annotated

from pydantic import BaseModel, Field, PositiveFloat, PositiveInt


class LoggingLevel(str, Enum):
    error = "ERROR"
    warning = "WARNING"
    info = "INFO"
    debug = "DEBUG"


class Placement(str, Enum):
    left = "left"
    right = "right"


class Layer(BaseModel):
    name: str
    zmin: float
    zmax: float
    index: Annotated[float, Field(ge=1.0)]


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
    placement: Placement = Placement.left
    yspan: PositiveFloat = 4.0
    zspan: PositiveFloat = 4.0
    y0: float = 0.0
    z0: float = 0.0
    port_resolution_y: PositiveInt = 30
    port_resolution_z: PositiveInt = 30


class Shape(BaseModel):
    cell_name: str = "cell"
    poly: List[Tuple[float, float]] = []
    xs_name: str = "default"


class ModeParams(BaseModel):
    logging_level: LoggingLevel = LoggingLevel.error
    eps_y: float = 0.0
    eps_z: float = 0.0
    std_y: PositiveFloat = 1.0
    std_z: PositiveFloat = 1.0
    max_iterations: PositiveInt = 1000
    min_iterations: PositiveInt = 10
    absolute_tolerance: PositiveFloat = 1.0e-4
    increment_x: PositiveFloat = 0.1
    get_increment_from_bpm: bool = True


class SettingsBPM(BaseModel):
    reference_index: Annotated[float, Field(ge=1.0)]
    wl: PositiveFloat
    resolution_x: PositiveInt
    resolution_y: PositiveInt
    resolution_z: PositiveInt
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
    scheme_parameter: Annotated[float, Field(ge=0.0, le=1.0)] = 0.5
    pml_strength: PositiveFloat = 5.0
    pml_thickness: PositiveFloat = 1.0
    dry_run: bool = False
    print_progress_percentage: List[Annotated[float, Field(gt=0.0, lt=100.0)]] = [
        1.0,
        10.0,
        20.0,
        30.0,
        40.0,
        50.0,
        100.0,
    ]
    index_slice_y: float = 0.0
    index_slice_z: float = 0.0
    field_slice_y: float = 0.0
    field_slice_z: float = 0.0
    mode_params: ModeParams
