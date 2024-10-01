Todo items:
1) Implement Crank-Nickelson scheme and compare speed and accuracy against RK4.
2) Implement Geometry class to keep track of polygons
3) Write to hdf5 instead of csv, and output more data to disk. Also read inputs from a file so you do not need to recompile every time.
4) Save also slices of the field at the beginning, middle, and end. For various directions. Add interpolation of the field so you can neatly save it at a given slice. And improve plotting.
5) Load mode field as initial profile, handle correct placement of the mode, and do the mode overlap at the output as well.
6) Compute mode fields inside the BPM itself.
7) Add unit tests.
8) Use Boost.MultiArray for multi-D arrays. This simplifies the array indexing, so it's less error-prone.
9) Add Python interface, swig, pybind11?
10) Add exception handling.