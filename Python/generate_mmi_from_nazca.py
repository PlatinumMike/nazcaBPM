#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:54:10 2024

MMI example

@author: mike
"""

import nazca as nd
from nazca.demofab import deep as ic
from nazca.geometries import box
from get_polygons_gds import extract_polygons_from_gds
import os
import json

WG_LAYER = (3, 0)


def slab(width: float = 4.0, length: float = 20.0):
    points = box(length=length, width=width)
    with nd.Cell("slab") as C:
        nd.Polygon(points=points, layer=WG_LAYER).put()
        nd.Pin("a0").put(0, width / 6, 180)
        nd.Pin("a1").put(0, -width / 6, 180)
        nd.Pin("b0").put(length, width / 6)
        nd.Pin("b1").put(length, -width / 6)
    return C


def get_mmi2x2(
    width_wg: float = 1.0,
    width_access: float = 1.5,
    width_slab: float = 5.0,
    length_slab: float = 10.0,
    taper_length: float = 2.0,
    strt_length: float = 5.0,
) -> nd.Cell:
    with nd.Cell("mmi") as C:
        mmi_slab = slab(width_slab, length_slab).put("org")
        for pin_name in ["a0", "a1", "b0", "b1"]:
            ic.ptaper(width1=width_access, width2=width_wg, length=taper_length).put(mmi_slab.pin[pin_name])
            port = ic.strt(length=strt_length, width=width_wg).put()
            port.raise_pins(["b0"], [pin_name])
    return C


def get_port_list(
    cell: nd.Cell,
    port_names: list,
    placement: str = "left",
    yspan: float = 2.0,
    zspan: float = 4.0,
    z0: float = 0.0,
    port_resolution_y: int = 20,
    port_resolution_z: int = 20,
) -> list:
    ports = []
    for port_name in port_names:
        port = {
            "name": port_name,
            "placement": placement,
            "yspan": yspan,
            "zspan": zspan,
            "y0": cell.pin[port_name].y,
            "z0": z0,
            "port_resolution_y": port_resolution_y,
            "port_resolution_z": port_resolution_z,
        }
        ports.append(port)
    return ports


width1 = 1.0
width2 = 1.8
taper_len = 10.0

# generate the gds
gds_file = "mmi_example.gds"
# note, the mmi length and width are arbitrary, this is not per realistic values.
mmi_cell = get_mmi2x2(width_wg=width1, width_access=width2, width_slab=8.0, length_slab=50.0, taper_length=4.0)
nd.export_gds(topcells=mmi_cell, filename=gds_file)


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


inports = get_port_list(
    cell=mmi_cell,
    port_names=["a0", "a1"],
    placement="left",
    yspan=2.0,
    zspan=height + 2,
    z0=0.0,
    port_resolution_y=resy,
    port_resolution_z=resz,
)

outports = get_port_list(
    cell=mmi_cell,
    port_names=["b0", "b1"],
    placement="right",
    yspan=2.0,
    zspan=height + 2,
    z0=0.0,
    port_resolution_y=resy,
    port_resolution_z=resz,
)


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
    "xmin": mmi_cell.pin["a0"].x + buffer,
    "xmax": mmi_cell.pin["b0"].x - buffer,
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
