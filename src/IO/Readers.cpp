//
// Created by mike on 10/1/24.
//

#include "Readers.h"
#include "../Geometry.h"
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
        read_json(inputFileName, root);

        //todo: split up struct into multiple pieces, not all these variables need to be known in every class.
        inputs.background_index = root.get<double>("background_index");
        inputs.reference_index = root.get<double>("reference_index");
        inputs.wl = root.get<double>("wl");
        inputs.resolution_x = root.get<int>("resolution_x");
        inputs.resolution_y = root.get<int>("resolution_y");
        inputs.resolution_z = root.get<int>("resolution_z");
        inputs.domain_len_x = root.get<double>("domain_len_x");
        inputs.domain_len_y = root.get<double>("domain_len_y");
        inputs.domain_len_z = root.get<double>("domain_len_z");
        inputs.pml_strength = root.get<double>("pml_strength");
        inputs.pml_thickness = root.get<double>("pml_thickness");
        inputs.scheme_parameter = root.get<double>("scheme_parameter");
        inputs.dry_run = root.get<bool>("dry_run");


        //compute derived quantities
        int numx = static_cast<int>(inputs.domain_len_x * inputs.resolution_x);
        int numy = static_cast<int>(inputs.domain_len_y * inputs.resolution_y);
        int numz = static_cast<int>(inputs.domain_len_z * inputs.resolution_z);


        inputs.numx = numx;
        inputs.numy = numy;
        inputs.numz = numz;
        inputs.num_slice = numx * numy;
        inputs.dx = inputs.domain_len_x / (numx - 1);
        inputs.dy = inputs.domain_len_y / (numy - 1);
        inputs.dz = inputs.domain_len_z / (numz - 1);


        inputs.k0 = 2 * std::numbers::pi / inputs.wl;
        inputs.beta_ref = inputs.k0 * inputs.reference_index;


        //read shapes, get max index

        double max_index = inputs.background_index;
        auto shapes = root.get_child("shapes");
        for (auto &shape: shapes) {
            auto actual_shape = shape.second; //get second element
            auto zmax = actual_shape.get<double>("zmax");
            auto zmin = actual_shape.get<double>("zmin");
            auto refractive_index = actual_shape.get<double>("refractive_index");
            if (refractive_index > max_index) {
                max_index = refractive_index;
            }
            auto poly_reader = actual_shape.get_child("poly");
            std::vector<Point> points;

            for (auto &array2: poly_reader) {
                double xposition, yposition;
                double *const elements[2] = {&xposition, &yposition};
                auto element = std::begin(elements);

                for (auto &i: array2.second) {
                    **element++ = i.second.get_value<double>();

                    if (element == std::end(elements)) break;
                }
                points.push_back({xposition, yposition});
            }
            Polygon poly(points);

            Shape new_shape(poly, zmin, zmax, refractive_index);
            inputs.shapes.push_back(new_shape);
        }
        inputs.max_index = max_index;
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
    assert(params.max_index >= params.background_index);
    assert(params.reference_index > 0);
    assert(params.wl > 0);
    assert(params.resolution_x > 0);
    assert(params.resolution_y > 0);
    assert(params.resolution_z > 0);
    assert(params.domain_len_x > 0);
    assert(params.domain_len_y > 0);
    assert(params.domain_len_z > 0);
    assert(params.pml_strength > 0);
    assert(params.pml_thickness > 0);
}
