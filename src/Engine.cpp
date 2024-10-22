//
// Created by mike on 9/29/24.
//

#include "Engine.h"
#include "IO/Readers.h"
#include "Solver.h"
#include "PML.h"
#include "ModeHandler.h"

#include <iostream>
#include <ostream>
#include <chrono>


Engine::Engine(const std::string &inputFileName) {
    std::cout << "Launching new engine now" << std::endl;
    //initialize
    _inputs = Readers::readJSON(inputFileName);

    std::cout << "Lx = " << _inputs.domain_len_x << ", Ly = " << _inputs.domain_len_y << ", Lz = " << _inputs.
            domain_len_z << std::endl;
    std::cout << "Nx = " << _inputs.numx << ", Ny = " << _inputs.numy << ", Nz = " << _inputs.numz << std::endl;
    std::cout << "dx = " << _inputs.dx << ", dy = " << _inputs.dy << ", dz = " << _inputs.dz << std::endl;
}

void Engine::run() const {
    const Parameters inputs = _inputs; //make a const copy, to avoid accidental modification of input variables.

    const double x_min = 0.0;
    const double x_max = inputs.domain_len_x;
    const double y_min = -0.5 * inputs.domain_len_y;
    const double y_max = 0.5 * inputs.domain_len_y;
    const double z_min = -0.5 * inputs.domain_len_z;
    const double z_max = 0.5 * inputs.domain_len_z;

    // define geometry
    const Geometry geometry(inputs.shapes, inputs.background_index);
    // set perfectly matched layers
    const PML pmly(inputs.pml_thickness, inputs.pml_strength, y_min, y_max); //adds a PML on both sides
    const PML pmlz(inputs.pml_thickness, inputs.pml_strength, z_min, z_max);
    // setup mode source
    const ModeHandler gaussianSource("Gaussian");


    Solver solver(geometry, pmly, pmlz, gaussianSource, x_min, x_max, y_min, y_max, z_min, z_max, inputs.numx,
                  inputs.numy, inputs.numz, inputs.scheme_parameter, inputs.k0, inputs.reference_index);
    solver.run();
    //todo: call function solver.save_data() or something to save/extract data after the solve is completed.
}
