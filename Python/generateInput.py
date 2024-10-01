"""
Generate input for BPM model.
"""

import json

# user inputs:
jsonFileName = "./inputs.json"

# convert to python dict
dataDict = {
    "background_index": 1.5,
    "core_index": 2.0,
    "reference_index": 1.6,
    "wl": 1.5,
    "resolution_x": 10,
    "resolution_y": 10,
    "domain_len_x": 4.0,
    "domain_len_y": 6.0,
    "domain_len_z": 10.0,
    "pml_strength": 5.0,
    "pml_thickness": 1.0,
}

# convert to JSON
jsonDict = json.dumps(dataDict, indent=2, sort_keys=True)

# write to disk
with open(jsonFileName, "w") as text_file:
    text_file.write(jsonDict)

# finished
print("json file created")
