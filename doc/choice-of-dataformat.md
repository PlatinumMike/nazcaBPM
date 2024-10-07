# Comparison of data formats

## CSV
It's easy to use, humanly readable, and you do not need any special libraries to read/write it.
The downside is the performance. E.g. storing an array of 1 billion doubles should only take 64 Gbit plus a tiny bit of overhead, but it takes quite a bit more with CSV because it's stored as ascii.
Also the read/write times are inferior to other methods.

## Raw binary dump
Actualy this is unformatted. And since you are not following any particular standard it is not guaranteed that the simulation results can even be read on another machine.
This has the highest performance, but it's terribly user unfriendly. Not recommended.

## HDF5
The hierachical data format sits somewhere in the middle. It is not as user friendly as CSV, but it has a far higher performance in both speed and file size. It is near the performance of unformatted IO, yet it uses a standardized format which ensures that it can still be read in the future. It is not humanly readable, but the tool `hdfview` can be used to inspect the files in a GUI. Also a python module exists for reading/writing HDF5 files. The [HDF group](https://www.hdfgroup.org/solutions/hdf5/) is here to stay it seems, so we can be sure that this file format will be supported for years to come. If not, we can always switch to something like [Zarr](https://zarr.dev/). And since the code is set-up in a modular way, switching to another kind of reader/writer should not be too much trouble. Also, you can continue to use this program in an isolated environment such as a docker long after official support for its libraries has ended.

## NetCDF
Netcdf is just HDF5 with extra stuff. Since we only need to store a few simple arrays plain HDF5 is good enough.
