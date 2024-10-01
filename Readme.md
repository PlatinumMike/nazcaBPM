# Beam propagation method
This simulation tool can be used to rapidly simulate structures, as long as they vary "slowly" in the propagation direction.
This includes: MMIs, directional couplers, adiabatic taper, star coupler, etc...
The simulation is performed in 3D, with perfectly matched layers (PMLs) in x and y direction.

Bends can be modelled as long as they are gradual enough (small angle). For 90 degree bends the propagation direction changes, so the slowly varying envelope approximation (SVEA) breaks down.
This can be tackled by using a conformal mapping, or changing to cylindrical coordinates. But currently the code is limited to simple Cartesian coordinates.

# Installation

You need `cmake` and `make`. Then just create a build directory:

`mkdir build; cd build`

Generate the Makefile:

`cmake .. -DCMAKE_BUILD_TYPE=Release`

And compile the program:

`make -j <ncores>`

with `<ncores>` the number of cores to compile with, e.g. 8.


# Geometry building

It is assumed that the geometry is given as a list of polygons in the (y,z) plane, and a height in x direction. The polygon is then simply extruded, so the side-wall angles are for now neglected.
To get the refractive index $n(x,y,z)$ at a given position, a temporary list is first constructed by looping over all shapes and checking if the current z position falls within the zmin and zmax of the shape.
Then, for that slice of the domain we calculate the index on the grid. For every point all shapes in the temporary list need to be evaluated to check if the current point is inside of the shape.
That can be done with the ray-casting algorithm. The list is traversed in reverse direction, so shapes that are added later take precedence.
This is important if there is an overlap. If the list is exhausted the (x,y) point does not fall within any of the shapes, so it is assigned the background index value.

# Misc

The code is not parallelized. But for the typical use case this is not an issue. We usually run large parameter sweeps of thousands of completely independent simulations.
These can run in parallel in separate processes. In other words, parallelize the sweep, not the individual simulation.
This should be highly efficient. Also, the Thomas algorithm does not parallelize well, so it is not trivial to get a speed-up inside a single simulation by using more cores.


# License
todo

Mike Machielsen 2024

This program comes with ABSOLUTELY NO WARRANTY.