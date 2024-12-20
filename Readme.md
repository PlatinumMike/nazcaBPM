# Beam propagation method
This simulation tool can be used to rapidly simulate structures, as long as they vary "slowly" in the propagation direction.
This includes: MMIs, directional couplers, adiabatic taper, star coupler, etc...
The simulation is performed in 3D, with perfectly matched layers (PMLs) in two transverse directions.

Bends can be modelled as long as they are gradual enough (small angle). For 90 degree bends the propagation direction changes, so the slowly varying envelope approximation (SVEA) breaks down.
This can be tackled by using a conformal mapping, or changing to cylindrical coordinates. But currently the code is limited to simple Cartesian coordinates.

# Prerequisites

To compile this program you need:
* C++ compiler (that supports c++20)
* cmake
* make (or ninja)
* HDF5 (On ubuntu this is easy, install `libhdf5-dev` with apt)
* Boost (install `libboost-all-dev`)
* python3 (optional)

Alternatively the program can run in a docker container, see below. 

## Conda
Another alternative way to install it is inside a conda environment (Warning: not yet tested). First install Anaconda or miniconda, then create a new environment, e.g. named `bpm`:

`conda create -n bpm`

Activate it:

`conda activate bpm`

Install dependencies:

`conda install -c conda-forge hdf5`

`conda install -c conda-forge gxx`

You can check with `which g++` to see if the compiler alias points towards your system's g++, or the one we just installed inside this environment. To compile the progam just continue the instructions in the next section. Conda forge has relatively fresh packages, so you're sure to have one of the latest compilers. This may be handy when your system's compiler is too old (need c++20).

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

# Docker
If you want to run the app in a docker, you do not need to install any prerequisites (except docker of course). To build the image from the docker file run the following in the same directory as where the `Dockerfile` is located.

`docker build -t nazca-bpm:0.2.0 .`

(do not forget the full-stop at the end.)
Then prepare a json input file as usual, and run the container with:

`docker run --rm -v ./:/mnt/workdir nazca-bpm:0.2.0 /mnt/workdir/Python/inputs.json`

The version can be specified after the colon, e.g. `0.2.0` in this case. The tag `nazca-bpm` is also abitrary, it is just a name. The `-v ./` mounts the current working directory to `/mnt/workdir` inside the container. So you can then specify where the json file is w.r.t. that. In this example your json input file should contain `"absolute_path_output": "/mnt/workdir"` so that in the docker it will write the data it produces into the mounted directory.

## Common issues:

Input file not found:
`Cannot find input file /mnt/workdir/Python/inputs.json`.
The docker mounts the current workding dir, and looks inside there for `Python/inputs.json`. So the json input file must be in your current working directory or deeper otherwise the docker cannot see it.

The container ran but there are no output files?
This is because the executable lives in `/app/build/` inside the container, but you only mount `/mnt/workdir/`. After the run stops any data of the container is removed. So make sure to specify the correct output dir in the json input file, so somewhere inside `/mnt/workdir/`.

# Geometry building

It is assumed that the geometry is given as a list of polygons in the (y,z) plane, and a height in x direction. The polygon is then simply extruded, so the side-wall angles are for now neglected.
To get the refractive index $n(x,y,z)$ at a given position, a temporary list is first constructed by looping over all shapes and checking if the current z position falls within the zmin and zmax of the shape.
Then, for that slice of the domain we calculate the index on the grid. For every point all shapes in the temporary list need to be evaluated to check if the current point is inside of the shape.
That can be done with the ray-casting algorithm. The list is traversed in reverse direction, so shapes that are added later take precedence.
This is important if there is an overlap. If the list is exhausted the (x,y) point does not fall within any of the shapes, so it is assigned the background index value.

Warning: Make sure that your geometry does not stop abruptly at the PML entrance, but fully extends into the layer. Otherwise you will still have reflections there. There is no performance penalty in letting the geometry extend beyond the modeling domain. In fact, it is recommended that you add tiny straight waveguides to the input/output ports of your device. This provides a good position for placing the mode source and monitors. Another caveat is that the model needs to check if a grid point is inside of a polygon, this can fail right at the edge of a polygon. So if your structure ends exactly at the domain boundary, it is unclear if that grid point will be inside or outside of your structure. If you then place your mode source in this exact slice which is outside of your your structure it cannot find any guided modes. So it is best to let the geometry extend a bit (1 nm is enough) beyond the domain boundaries in propagation direction.

# Testing

The software ships with some basic test cases that you can run. These are opt-in, to enable them first download the [doctest](https://github.com/doctest/doctest) header ``doctest.h`` and place it inside the test directory. Secondly, uncomment the `set(ENABLE_TEST 1)` in `CMakeLists.txt`, go to the build directory and run `cmake` followed by `make`. Lastly, the tests can be run with ``./testSuite.x``. Tested so far on Ubuntu 24.04 with gcc 13.2.0 and boost 1.83.0.

# Misc

The code is not parallelized. But for the typical use case this is not an issue. We usually run large parameter sweeps of thousands of completely independent simulations.
These can run in parallel in separate processes. In other words, parallelize the sweep, not the individual simulation.
This should be highly efficient. Also, the Thomas algorithm does not parallelize well, so it is not trivial to get a speed-up inside a single simulation by using more cores.


# License
Copyright (c) 2024 Mike Machielsen

See `LICENSE` file for the details.

