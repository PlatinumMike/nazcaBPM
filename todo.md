Todo items:
- Save also slices of the field at the beginning, middle, and end. For various directions. Add interpolation of the field so you can neatly save it at a given slice. And improve plotting.
- Load mode field as initial profile, handle correct placement of the mode, and do the mode overlap at the output as well.
 The mode solver uses the monitor grid, which is not the exact same as the BPM main grid. This results in some small oscillations in x direction before the mode settles on the BPM grid.
 Can we fix this?
- Compute higher order modes inside the BPM.
- Add more unit tests.
- Add Python interface, swig, pybind11? Or dockerize?
- Add exception handling.
- Implement a semi-vectorial version.
- Implement conformal mapping, to handle bends? If you cut up the bend into a series of circular arc segments you can model any kind of spline, including cobras and euler bends.
 It's not scalable to have to cut up your polygons into many small pieces. Then checking if a point is inside a polygon becomes more work because there are more polygons.
 A better way would be to add another entry in the input file where you specify the curvature as a function of propagation distance x. This will then apply to any polygon.
 During evaluation of the index it will first compute the normal index and then modify it using the local curvature. The curvature can be any function of x, so it's probably
 easiest to just ask the user to give a range of x points with their corresponding curvature value. In the code we then use 1D linear interpolation to get the curvature anywhere.
- Add folder with documentation on the derivation of the method, and its pros and cons over other simulation techniques.
- Once the code is more mature, add more docs on how to use it. Use Sphinx or mkdocs.
- Make it optional to export the refractive index and field in slices. This is useful initially when setting up a simulation to verify that it models what you think it does. But it does not make sense to export this for all 10^5 simulations in the sweep. So for the other 99999 simulations just export the mode overlaps at the output. If you are worried about flying blind and modeling the wrong thing, just rerun the final optimized scenario with output of the fields to check.
- Profile the code to see where it can be sped up.
- Translate to Rust? Replace boost by `ndarray` crate. For json use `serde`. We can use `hdf5`, or use `ndarray` to .npy files (numpy). If you make a proper python interface we can even avoid writing to disk entirely, the results just stay in memory. If you eliminate hdf5, you also do not need any external libs, so it will be easier to wrap it as a python package.
- Add some simple script sim.plot() that shows the geometry and PML so that a user knows the setup is correct. Similar to what MEEP does.

