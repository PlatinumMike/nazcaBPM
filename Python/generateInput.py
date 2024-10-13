"""
Generate input for BPM model.
"""

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
        (0.1, width1 / 2),
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


# convert to python dict
dataDict = {
    "background_index": 1.5,
    "reference_index": 1.6,
    "wl": 1.5,
    "resolution_x": 10,
    "resolution_y": 10,
    "resolution_z": 10,
    "domain_len_x": 4.0,
    "domain_len_y": 6.0,
    "domain_len_z": taper_length + 30.0,
    "pml_strength": 5.0,
    "pml_thickness": 1.0,
    "shapes": shapes,
    "scheme_parameter": 0.5,
}

# convert to JSON
jsonDict = json.dumps(dataDict, indent=2, sort_keys=True)

# write to disk
with open(jsonFileName, "w") as text_file:
    text_file.write(jsonDict)

# finished
print("json file created")
