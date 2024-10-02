//
// Created by mike on 10/1/24.
//

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Geometry.h"

struct Parameters {

    double background_index; //background refractive index
    double max_index; //maximum refractive index in the entire domain
    double wl; //vacuum wavelength (um)
    int resolution_x; //number of gridpoints per unit length in x direction
    int resolution_y; //number of gridpoints per unit length in y direction
    int numx; //number of gridpoints
    int numy;
    int numz;
    int num_slice; // number of gridpoints in one xy slice
    double domain_len_x; //length of the domain in x direction
    double domain_len_y; //length of the domain in y direction
    double domain_len_z; //length of the domain in z direction
    double reference_index;
    double k0; //vacuum wavenumber
    double beta_ref;

    double pml_strength; //max conductivity inside the PML
    double pml_thickness; //thickness of the PML (um)

    std::vector<Shape> shapes;

    //todo: add method to print all of the loaded values inside of the code, so a user can confirm they are loaded correctly.
};

#endif //PARAMETERS_H
