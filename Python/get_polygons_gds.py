#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 20:09:04 2024

test reading a gds.
Just make sure to flatten the layout, this way we are sure there is only one cell to load.
Made with chatgpt.

Note, this only seems to work for binary GDS2, not ascii data.

@author: mike
"""


import gdstk

# TODO: layer number needs to be an integer, but what if we pass in a string or tuple?


# Function to extract polygons from a GDS file
def extract_polygons_from_gds(file_path, layer_number: int):
    # Read the GDS file
    library = gdstk.read_gds(file_path)
    polygons = []

    # Iterate over all cells in the GDS library
    for cell in library.cells:
        # Iterate over all elements in the cell
        for element in cell.polygons:
            # Append the polygon coordinates to the list
            if element.layer == layer_number:
                polygons.append(element.points)

    return polygons


if __name__ == "__main__":
    # Example usage
    gds_file = "example.gds"  # Replace with the path to your GDS file
    polygons = extract_polygons_from_gds(gds_file, 1)

    # Print the extracted polygons
    print(f"Extracted {len(polygons)} polygons:")
    for i, polygon in enumerate(polygons, 1):
        print(f"Polygon {i}: {polygon}")
