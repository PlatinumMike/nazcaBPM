# Beam propagation method
This simulation tool can be used to rapidly simulate structures, as long as they vary "slowly" in the propagation direction.
This includes: MMIs, directional couplers, adiabatic taper, star coupler, etc...
The simulation is performed in 3D, with perfectly matched layers (PMLs) in x and y direction.

Bends can be modelled as long as they are gradual enough (small angle). For 90 degree bends the propagation direction changes, so the slowly varying envelope approximation (SVEA) breaks down.
This can be tackled by using a conformal mapping, or changing to cylindrical coordinates. But currently the code is limited to simple Cartesian coordinates.

# Prerequisites

To compile this program you need:
* C++ compiler
* cmake
* make (or ninja)
* HDF5 (On ubuntu this is easy, install `libhdf5-dev`)
* Boost (install `libboost-all-dev`)
* python3 (optional)

# Installation

You need `cmake` and `make`. Then just create a build directory:

`mkdir build; cd build`

Generate the Makefile:

`cmake .. -DCMAKE_BUILD_TYPE=Release`

And compile the program:

`make -j <ncores>`

with `<ncores>` the number of cores to compile with, e.g. 8.

It may happen that cmake cannot find the required libraries. If you install them with `apt` they will be added to the include path automatically. 
If you just download a zip file, `cmake` will of course have no idea where to look, unless you manually add this to the include path. 

# Running the program
The executable in inside the build directory. It takes one command line argument, the name of the input file.
This input file can be generated by running generateInput.py. Example inputs are included in the Inputs directory.

`./nazcaBPM.x inputs.json`

# Geometry building

It is assumed that the geometry is given as a list of polygons in the (y,z) plane, and a height in x direction. The polygon is then simply extruded, so the side-wall angles are for now neglected.
To get the refractive index $n(x,y,z)$ at a given position, a temporary list is first constructed by looping over all shapes and checking if the current z position falls within the zmin and zmax of the shape.
Then, for that slice of the domain we calculate the index on the grid. For every point all shapes in the temporary list need to be evaluated to check if the current point is inside of the shape.
That can be done with the ray-casting algorithm. The list is traversed in reverse direction, so shapes that are added later take precedence.
This is important if there is an overlap. If the list is exhausted the (x,y) point does not fall within any of the shapes, so it is assigned the background index value.

# Testing

The software ships with some basic test cases that you can run. These are opt-in, to enable them first download the [doctest](https://github.com/doctest/doctest) header ``doctest.h`` and place it inside the test directory. Secondly, uncomment the `set(ENABLE_TEST 1)` in `CMakeLists.txt`, go to the build directory and run `cmake` followed by `make`. Lastly, the tests can be run with ``./testSuite.x``. Tested so far on Ubuntu 24.04 with gcc 13.2.0 and boost 1.83.0.

# Misc

The code is not parallelized. But for the typical use case this is not an issue. We usually run large parameter sweeps of thousands of completely independent simulations.
These can run in parallel in separate processes. In other words, parallelize the sweep, not the individual simulation.
This should be highly efficient. Also, the Thomas algorithm does not parallelize well, so it is not trivial to get a speed-up inside a single simulation by using more cores.


# License
Copyright (c) 2024 Mike Machielsen

See `LICENSE` file for the details.

