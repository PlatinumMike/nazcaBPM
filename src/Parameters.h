//
// Created by mike on 10/1/24.
//

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <filesystem>
#include "Geometry.h"
#include "Port.h"

struct Parameters {
    double wl; //vacuum wavelength (um)
    int resolution_x; //number of gridpoints per unit length in x direction
    int resolution_y; //number of gridpoints per unit length in y direction
    int resolution_z; //number of gridpoints per unit length in z direction
    int numx; //number of gridpoints
    int numy;
    int numz;
    int num_slice; // number of gridpoints in one xy slice
    double domain_len_x; //length of the domain in x direction
    double domain_len_y; //length of the domain in y direction
    double domain_len_z; //length of the domain in z direction
    double xmin; // x position of the domain, lower bound
    double xmax; // x position of the domain, upper bound
    double ymin; // y position of the domain, lower bound
    double ymax; // y position of the domain, upper bound
    double zmin; // z position of the domain, lower bound
    double zmax; // z position of the domain, upper bound
    double dx; //grid spacing, assuming a uniform mesh
    double dy;
    double dz;
    double reference_index;
    double k0; //vacuum wavenumber
    double beta_ref;

    double scheme_parameter; //value between 0 and 1. For the classical Crank-Nicolson use 0.5.

    double pml_strength; //max conductivity inside the PML
    double pml_thickness; //thickness of the PML (um)

    // If dry_run is true, just load the parameters, evaluate refractive index on the grid and exit the program.
    // This is useful before launching a big simulation. You want to know if you are modeling the correct thing before launching into a long run.
    bool dry_run;

    std::filesystem::path absolute_path_output; //absolute path of where the generated files are to be stored.


    std::vector<Shape> shapes;
    std::vector<Port> input_ports;
    std::vector<Port> output_ports;

    std::unordered_map<std::string, XS> xs_map;

    //todo: add method to print all of the loaded values inside of the code, so a user can confirm they are loaded correctly.
};

#endif //PARAMETERS_H
