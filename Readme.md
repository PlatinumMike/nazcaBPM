# Beam propagation method
This 3D simulation tool can be used to rapidly simulate structures, as long as they vary "slowly" in the propagation direction.
This includes: MMIs, directional couplers, adiabatic taper, star coupler, etc...

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
todo

# License
todo

Mike Machielsen 2024

This program comes with ABSOLUTELY NO WARRANTY.