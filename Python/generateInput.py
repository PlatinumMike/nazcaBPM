"""
Generate input for BPM model.
"""

import os
import json

# user inputs:
jsonFileName = "./inputs.json"

# simple taper geometry

width1 = 1.000
width2 = 3.000
height = 0.350
taper_length = 10
zmin = -height / 2
zmax = height / 2

resx = 10
resy = 30
resz = resy

# this polygon could be exported by nazca for instance.
# Adding extra bit of straight waveguide, otherwise the polygon ends right on the domain border, which means there is a jump in index at the start.
taper = {
    "cell_name": "taper",
    "zmin": zmin,
    "zmax": zmax,
    "refractive_index": 2.0,
    "poly": [
        (-0.1, -width1 / 2),
        (10, -width1 / 2),
        (10 + taper_length, -width2 / 2),
        (20 + taper_length, -width2 / 2),
        (20 + taper_length, width2 / 2),
        (10 + taper_length, width2 / 2),
        (10, width1 / 2),
        (-0.1, width1 / 2),
    ],
}
# adding a straight waveguide
strt = {
    "cell_name": "strt",
    "zmin": zmin,
    "zmax": zmax,
    "refractive_index": 2.0,
    "poly": [
        (20 + taper_length, -width2 / 2),
        (30.1 + taper_length, -width2 / 2),
        (30.1 + taper_length, width2 / 2),
        (20 + taper_length, width2 / 2),
    ],
}
shapes = [taper, strt]

# ports have their own grid, and thus their own resolution.
# no point in using a higher resolution than that of the main simulation because the field is interpolated on the BPM simulation grid anyway.
port_a0 = {
    "name": "a0",
    "placement": "left",
    "yspan": width1 + 2,
    "zspan": height + 2,
    "y0": 0.0,
    "z0": 0.0,
    "port_resolution_y": resy,
    "port_resolution_z": resz,
}
port_b0 = {
    "name": "b0",
    "placement": "right",
    "yspan": width2 + 2,
    "zspan": height + 2,
    "y0": 0.0,
    "z0": 0.0,
    "port_resolution_y": resy,
    "port_resolution_z": resz,
}

inports = [port_a0]
outports = [port_b0]


# convert to python dict
dataDict = {
    "background_index": 1.5,
    "reference_index": 1.65,
    "wl": 1.5,
    "resolution_x": resx,
    "resolution_y": resy,
    "resolution_z": resz,
    "xmin": 0.0,
    "xmax": taper_length + 30.0,
    "ymin": -3.0,
    "ymax": 3.0,
    "zmin": -2.0,
    "zmax": 2.0,
    "pml_strength": 5.0,
    "pml_thickness": 1.0,
    "shapes": shapes,
    "scheme_parameter": 0.5,
    "dry_run": False,
    "input_ports": inports,
    "output_ports": outports,
}

# convert to JSON
jsonDict = json.dumps(dataDict, indent=2, sort_keys=True)

# write to disk
fname = os.path.abspath(jsonFileName)
with open(fname, "w") as text_file:
    text_file.write(jsonDict)

# finished
print(f"json file {fname} created")
