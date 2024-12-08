#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:54:10 2024

Get polygons from GDS file, then make an input file for the BPM simulation.
In principle you could extract the polygons from a nazca cell directly.
But often you want to export a gds file anyway so you can easily examine the layout.
It is then convient to just import that gds file.

Here we show a simple taper example.

@author: mike
"""

import nazca as nd
from nazca.demofab import deep as ic
from get_polygons_gds import extract_polygons_from_gds
import os
import json


def get_taper(
    width1: float = 1.0, width2: float = 3.0, taper_len: float = 20.0, strt_len: float = 5.0
) -> nd.Cell:
    with nd.Cell("my-taper") as C:
        start = ic.strt(length=strt_len, width=width1).put()
        ic.ptaper(length=taper_len, width1=width1, width2=width2).put()
        end = ic.strt(length=strt_len, width=width2).put()

        # raise pins
        start.raise_pins(["a0"])
        end.raise_pins(["b0"])
    return C


width1 = 1.0
width2 = 3.0
taper_len = 10.0

# generate the gds
gds_file = "taper_example.gds"
taper_cell = get_taper(width1=width1, width2=width2, taper_len=taper_len)
nd.export_gds(topcells=taper_cell, filename=gds_file)


# load polygons
polygons = extract_polygons_from_gds(gds_file, 3)


this_dir = os.path.dirname(__file__)
root = os.path.dirname(this_dir)
working_dir = os.path.join(root, "simulations")

# other inputs
height = 0.350
zmin = -height / 2
zmax = height / 2

resx = 10
resy = 30
resz = resy

shapes = []
for i, polygon in enumerate(polygons):
    shape = {
        "cell_name": f"cell{i}",
        "zmin": zmin,
        "zmax": zmax,
        "refractive_index": 2.0,
        "poly": polygon.tolist(),
    }
    shapes.append(shape)


port_a0 = {
    "name": "a0",
    "placement": "left",
    "yspan": width1 + 2,
    "zspan": height + 2,
    "y0": taper_cell.pin["a0"].y,
    "z0": 0.0,
    "port_resolution_y": resy,
    "port_resolution_z": resz,
}
port_b0 = {
    "name": "b0",
    "placement": "right",
    "yspan": width2 + 2,
    "zspan": height + 2,
    "y0": taper_cell.pin["b0"].y,
    "z0": 0.0,
    "port_resolution_y": resy,
    "port_resolution_z": resz,
}

inports = [port_a0]
outports = [port_b0]

# move in the starting and ending boundaries a bit to ensure the structure extends all the way through the xmin, xmax
buffer = 1.0

# convert to python dict
dataDict = {
    "background_index": 1.5,
    "reference_index": 1.65,
    "wl": 1.55,
    "resolution_x": resx,
    "resolution_y": resy,
    "resolution_z": resz,
    "xmin": taper_cell.pin["a0"].x + buffer,
    "xmax": taper_cell.pin["b0"].x - buffer,
    "ymin": -8.0,
    "ymax": 8.0,
    "zmin": -4.0,
    "zmax": 4.0,
    "pml_strength": 5.0,
    "pml_thickness": 1.0,
    "shapes": shapes,
    "scheme_parameter": 0.5,
    "dry_run": False,
    "input_ports": inports,
    "output_ports": outports,
    "absolute_path_output": working_dir,
}

# convert to JSON
jsonDict = json.dumps(dataDict, indent=2, sort_keys=True)

jsonFileName = "./inputs.json"

# write to disk
fname = os.path.abspath(jsonFileName)
with open(fname, "w") as text_file:
    text_file.write(jsonDict)

# finished
print(f"json file {fname} created")
