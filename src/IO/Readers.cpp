//
// Created by mike on 10/1/24.
//

#include "Readers.h"
#include <iostream>
#include <filesystem>
#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"
#include <cassert>


Parameters Readers::readJSON(const std::string &inputFileName) {
    Parameters inputs{};
    std::cout << "Reading input from file " << inputFileName << std::endl;
    if (std::filesystem::exists(inputFileName)) {
        boost::property_tree::ptree root;
        boost::property_tree::read_json(inputFileName, root);

        //todo: split up struct into multiple pieces, not all these variables need to be known in every class.
        inputs.background_index = root.get<double>("background_index");
        inputs.core_index = root.get<double>("core_index");
        inputs.reference_index = root.get<double>("reference_index");
        inputs.wl = root.get<double>("wl");
        inputs.resolution_x = root.get<int>("resolution_x");
        inputs.resolution_y = root.get<int>("resolution_y");
        inputs.domain_len_x = root.get<double>("domain_len_x");
        inputs.domain_len_y = root.get<double>("domain_len_y");
        inputs.domain_len_z = root.get<double>("domain_len_z");
        inputs.pml_strength = root.get<double>("pml_strength");
        inputs.pml_thickness = root.get<double>("pml_thickness");


        //compute derived quantities
        int numx = static_cast<int>(inputs.domain_len_x * inputs.resolution_x);
        int numy = static_cast<int>(inputs.domain_len_y * inputs.resolution_y);


        inputs.numx = numx;
        inputs.numy = numy;
        inputs.num_slice = numx * numy;
        inputs.k0 = 2 * std::numbers::pi / inputs.wl;
        inputs.beta_ref = inputs.k0 * inputs.reference_index;
    } else {
        std::cout << "Cannot find input file " << inputFileName << std::endl;
        exit(EXIT_FAILURE);
    }

    testParameters(inputs);

    return inputs;
}

void Readers::testParameters(const Parameters &params) {
    //todo: asserts are ignored in release mode, find a way to enable them.
    assert(params.background_index > 0);
    assert(params.core_index > 0);
    assert(params.reference_index > 0);
    assert(params.wl > 0);
    assert(params.resolution_x > 0);
    assert(params.resolution_y > 0);
    assert(params.domain_len_x > 0);
    assert(params.domain_len_y > 0);
    assert(params.domain_len_z > 0);
    assert(params.pml_strength > 0);
    assert(params.pml_thickness > 0);
}
