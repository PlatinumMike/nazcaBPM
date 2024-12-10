//
// Created by mike on 10/1/24.
//

#include "Readers.h"
#include "../Geometry.h"
#include "../Port.h"
#include "../AuxiliaryFunctions.h"
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
    auto abs_path = root.get<std::string>("absolute_path_output");
    inputs.absolute_path_output = abs_path;


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


    //read cross-sections
    inputs.xs_map = get_xs_map(root);
    //read shapes
    inputs.shapes = get_shapes(root);
    // read in and output ports
    inputs.input_ports = get_ports(root, "input_ports", inputs.xmin, inputs.xmax, inputs.ymin, inputs.ymax, inputs.zmin,
                                   inputs.zmax);
    inputs.output_ports = get_ports(root, "output_ports", inputs.xmin, inputs.xmax, inputs.ymin, inputs.ymax,
                                    inputs.zmin, inputs.zmax);

    testParameters(inputs);

    return inputs;
}

void Readers::testParameters(const Parameters &params) {
    //todo: asserts are ignored in release mode, find a way to enable them.
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
    // the map needs to have a default cross-section in case the sampled point lies outside of a polygon in your GDS.
    assert(params.xs_map.contains("default"));
}

std::vector<Port> Readers::get_ports(boost::property_tree::ptree root, const std::string &portnames, const double xmin,
                                     const double xmax, double ymin, double ymax, double zmin, double zmax) {
    std::vector<Port> portVec;
    auto ports = root.get_child(portnames);
    for (auto &port: ports) {
        auto actual_port = port.second; //get second element
        const auto y0 = actual_port.get<double>("y0");
        const auto yspan = actual_port.get<double>("yspan");
        const auto z0 = actual_port.get<double>("z0");
        const auto zspan = actual_port.get<double>("zspan");
        const auto resolution_y = actual_port.get<double>("port_resolution_y");
        const auto resolution_z = actual_port.get<double>("port_resolution_z");
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
        double ymin_port = y0 - 0.5 * yspan;
        double ymax_port = y0 + 0.5 * yspan;
        double zmin_port = z0 - 0.5 * zspan;
        double zmax_port = z0 + 0.5 * zspan;

        //check if port grid does not extend beyond the total simulation grid
        ymin_port = std::max(ymin_port, ymin);
        ymax_port = std::min(ymax_port, ymax);
        zmin_port = std::max(zmin_port, zmin);
        zmax_port = std::min(zmax_port, zmax);
        //recompute yspan, zspan in case those have changed.
        const double yspan_port = ymax_port - ymin_port;
        const double zspan_port = zmax_port - zmin_port;

        int numy_port = static_cast<int>(yspan_port * resolution_y);
        int numz_port = static_cast<int>(zspan_port * resolution_z);
        auto ygrid = AuxiliaryFunctions::linspace(ymin_port, ymax_port, numy_port);
        auto zgrid = AuxiliaryFunctions::linspace(zmin_port, zmax_port, numz_port);
        Port new_port(name, placement, x0, ygrid, zgrid);
        portVec.push_back(new_port);
    }

    return portVec;
}

std::vector<Shape> Readers::get_shapes(boost::property_tree::ptree root) {
    std::vector<Shape> shapesVec;
    auto shapes = root.get_child("shapes");
    for (auto &shape: shapes) {
        auto actual_shape = shape.second; //get second element
        const auto xs_name = actual_shape.get<std::string>("xs_name");
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

        Shape new_shape(points, xs_name);
        shapesVec.push_back(new_shape);
    }
    return shapesVec;
}


std::unordered_map<std::string, XS> Readers::get_xs_map(boost::property_tree::ptree root) {
    std::unordered_map<std::string, XS> xs_map;
    auto cross_sections = root.get_child("cross_sections");
    for (auto &xs: cross_sections) {
        auto actual_xs = xs.second; //get second element
        const auto xs_name = actual_xs.get<std::string>("xs_name");
        const auto background_index = actual_xs.get<double>("background_index");
        auto list_reader = actual_xs.get_child("layer_list");

        XS new_xs(background_index);

        for (auto &layer: list_reader) {
            auto actual_layer = layer.second;
            const auto index = actual_layer.get<double>("index");
            const auto zmin = actual_layer.get<double>("zmin");
            const auto zmax = actual_layer.get<double>("zmax");
            const Layer new_layer(zmin, zmax, index);
            new_xs.append_layer(new_layer);
        }
        xs_map.insert({xs_name, new_xs});
    }
    return xs_map;
}
