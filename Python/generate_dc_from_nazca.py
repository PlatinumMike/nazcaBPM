#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:54:10 2024

Directional coupler example.
Note, the BPM works best for propagation close to the x-axis, so for small angles.
Because of these bends the field calculated by BPM will not be perfectly accurate. The tighter the bends, the worse.
It is advised to only use BPM as a tool to rapidly explore the parameter space, then use FDTD to do the final optimization.

@author: mike
"""

import nazca as nd
from nazca.demofab import deep as ic
from get_polygons_gds import extract_polygons_from_gds
import os
import json
from xs_information import sin_strip


def get_dc(
    width_wg: float = 1.0, offset: float = 10.0, gap: float = 0.3, length: float = 50.0, strt_len: float = 5.0
):
    with nd.Cell("directional-coupler") as C:
        ina0 = ic.sbend(offset=-offset, width=width_wg, length1=strt_len, length2=length).put(
            0, offset + 0.5 * (width_wg + gap)
        )
        outb0 = ic.sbend(offset=offset, width=width_wg, length2=strt_len).put()

        ina1 = ic.sbend(offset=offset, width=width_wg, length1=strt_len, length2=length).put(
            0, -offset - 0.5 * (width_wg + gap)
        )
        outb1 = ic.sbend(offset=-offset, width=width_wg, length2=strt_len).put()

        ina0.raise_pins(["a0"], ["a0"])
        ina1.raise_pins(["a0"], ["a1"])
        outb0.raise_pins(["b0"], ["b0"])
        outb1.raise_pins(["b0"], ["b1"])
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


width = 1.0
gap = 0.3

# generate the gds
gds_file = "dc_example.gds"
# note, the widths and gap are chosen arbitrarily, not optimized for 3dB coupling.
dc_cell = get_dc(width_wg=width, offset=8, gap=gap)
nd.export_gds(topcells=dc_cell, filename=gds_file)


# load polygons
polygons = extract_polygons_from_gds(gds_file, 3)


this_dir = os.path.dirname(__file__)
root = os.path.dirname(this_dir)
working_dir = os.path.join(root, "simulations")

# other inputs
height = 0.350

resx = 10
resy = 30
resz = resy

shapes = []
for i, polygon in enumerate(polygons):
    shape = {"cell_name": f"cell{i}", "poly": polygon.tolist(), "xs_name": "sin-strip"}
    shapes.append(shape)

xs_core, xs_default = sin_strip(height=height)

cross_sections = [xs_core.xs2str(), xs_default.xs2str()]


inports = get_port_list(
    cell=dc_cell,
    port_names=["a0", "a1"],
    placement="left",
    yspan=2.0,
    zspan=height + 2,
    z0=0.0,
    port_resolution_y=resy,
    port_resolution_z=resz,
)

outports = get_port_list(
    cell=dc_cell,
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
    "reference_index": 1.65,
    "wl": 1.55,
    "resolution_x": resx,
    "resolution_y": resy,
    "resolution_z": resz,
    "xmin": dc_cell.pin["a0"].x + buffer,
    "xmax": dc_cell.pin["b0"].x - buffer,
    "ymin": -12.0,
    "ymax": 12.0,
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
    "cross_sections": cross_sections,
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
