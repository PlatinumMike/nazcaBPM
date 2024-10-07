# Getting started

Read this before you start simulating.

## What is BPM
The beam propagation method (BPM) is a way of solving Maxwell's equations in the frequency domain. It makes the slowly varying envelope approximation (SVEA), which means the $E, H$ fields are written as the product of an envelope function and a rapidly varying phase part. The envelope varies only slowly in the propagation direction, compared to the wavelength scale. This is particularly well suited for modelling very long adiabatic structures, such as tapers, MMIs and directional couplers.

### Pros
1. It's very fast, typically taking just a few seconds for one simulation.
2. Finite difference BPM is relatively simple to implement.

### Cons
1. It omits reflections (unless you use bi-directional BPM)
2. The paraxial approximation means it only works when the light propagates close to the optical axis (although you could use wide-angle BPM to mitigate this)
3. The structures in the model must be slowly varying in the propagation direction.

### Conclusion

For a rigorous simulation you have to use EME or FDTD. But that does not mean BPM is useless, far from it. BPM can be used to rapidly explore the parameter space (when the SVEA is justified), leaving the number crunching of EME/FDTD for just the final simulations.

## Installation
See Readme.md

## First simulation
We need to prepare an input file for the simulation to work with. Run `generateInput.py` to generate one, or use the example in the `inputs` directory. Copy this to the build directory, then run `./nazcaBPM example_inputs.json`. Give it a few seconds to complete. The output is written as HDF5 files. This format is highly suitable for reading/writing large (multi-dimensional) arrays. 

The next step is to prepare the python environment. It is highly recommended to create a new environment to work in, so you don't pollute your system python with all kinds of packages. To do this run `python -m venv my_new_env` and activate it with `source my_new_env/bin/activate`. Then navigate to the Python directory and install the packages with `pip install -r requirements.txt`.

You can now plot the results of your first simulation with `plot_sim.py`.