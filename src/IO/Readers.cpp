//
// Created by mike on 10/1/24.
//

#include "Readers.h"
#include "../Geometry.h"
#include "../Port.h"
#include <iostream>
#include <filesystem>
#include <format>
#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"
#include <cassert>


Parameters Readers::readJSON(const std::string &inputFileName) {
    Parameters inputs{};
    std::cout << "Reading input from file " << inputFileName << std::endl;

    if (!std::filesystem::exists(inputFileName)) {
        std::cout << "Cannot find input file " << inputFileName << std::endl;
        exit(EXIT_FAILURE);
    }

    boost::property_tree::ptree root;
    read_json(inputFileName, root);

    //todo: split up struct into multiple pieces, not all these variables need to be known in every class.
    inputs.background_index = root.get<double>("background_index");
    inputs.reference_index = root.get<double>("reference_index");
    inputs.wl = root.get<double>("wl");
    inputs.resolution_x = root.get<int>("resolution_x");
    inputs.resolution_y = root.get<int>("resolution_y");
    inputs.resolution_z = root.get<int>("resolution_z");
    inputs.xmin = root.get<double>("xmin");
    inputs.xmax = root.get<double>("xmax");
    inputs.ymin = root.get<double>("ymin");
    inputs.ymax = root.get<double>("ymax");
    inputs.zmin = root.get<double>("zmin");
    inputs.zmax = root.get<double>("zmax");
    inputs.pml_strength = root.get<double>("pml_strength");
    inputs.pml_thickness = root.get<double>("pml_thickness");
    inputs.scheme_parameter = root.get<double>("scheme_parameter");
    inputs.dry_run = root.get<bool>("dry_run");


    inputs.domain_len_x = inputs.xmax - inputs.xmin;
    inputs.domain_len_y = inputs.ymax - inputs.ymin;
    inputs.domain_len_z = inputs.zmax - inputs.zmin;

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


    //read shapes
    inputs.shapes = get_shapes(root);
    //get max index
    inputs.max_index = get_max_index(inputs.shapes, inputs.background_index);
    // read in and output ports
    inputs.input_ports = get_ports(root, "input_ports", inputs.xmin, inputs.xmax);
    inputs.output_ports = get_ports(root, "output_ports", inputs.xmin, inputs.xmax);

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

std::vector<Port> Readers::get_ports(boost::property_tree::ptree root, const std::string &portnames, const double xmin,
                                     const double xmax) {
    std::vector<Port> portVec;
    auto ports = root.get_child(portnames);
    for (auto &port: ports) {
        auto actual_port = port.second; //get second element
        const auto y0 = actual_port.get<double>("y0");
        const auto yspan = actual_port.get<double>("yspan");
        const auto z0 = actual_port.get<double>("z0");
        const auto zspan = actual_port.get<double>("zspan");
        const auto name = actual_port.get<std::string>("name");
        const auto placement = actual_port.get<std::string>("placement");
        double x0 = 0.0;
        if (placement == "left") {
            x0 = xmin;
        } else if (placement == "right") {
            x0 = xmax;
        } else {
            std::cout << std::format("Port {} placed incorrectly! Placement needs to be 'right' or 'left'.", name)
                    << std::endl;
            exit(EXIT_FAILURE);
        }
        Port new_port(name, placement, x0, y0, z0, yspan, zspan);
        portVec.push_back(new_port);
    }

    return portVec;
}

std::vector<Shape> Readers::get_shapes(boost::property_tree::ptree root) {
    std::vector<Shape> shapesVec;
    auto shapes = root.get_child("shapes");
    for (auto &shape: shapes) {
        auto actual_shape = shape.second; //get second element
        auto zmax = actual_shape.get<double>("zmax");
        auto zmin = actual_shape.get<double>("zmin");
        auto refractive_index = actual_shape.get<double>("refractive_index");
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
        const Polygon poly(points);

        Shape new_shape(poly, zmin, zmax, refractive_index);
        shapesVec.push_back(new_shape);
    }
    return shapesVec;
}

double Readers::get_max_index(const std::vector<Shape> &shapes, const double background_index) {
    double max_index = background_index;
    for (const auto &shape: shapes) {
        max_index = std::max(max_index, shape.refractive_index);
    }
    return max_index;
}
