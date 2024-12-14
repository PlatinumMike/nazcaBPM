//
// Created by mike on 10/1/24.
//

#ifndef READERS_H
#define READERS_H

#include "../Parameters.h"
#include <string>
#include "boost/property_tree/ptree.hpp"


class Readers {
public:
    static Parameters readJSON(const std::string &inputFileName);

private:
    /**
     * Do some basic tests on the inputs before launching a simulation.
     * This avoids wasting large amounts of compute time on invalid inputs.
     * @param params input parameters
     */
    static void testParameters(const Parameters &params);

    static std::vector<Port> get_ports(boost::property_tree::ptree root, const std::string& portnames, double xmin,
                                       double xmax, double ymin, double ymax, double zmin, double zmax);

    static std::vector<Shape> get_shapes(boost::property_tree::ptree root);

    static std::unordered_map<std::string, XS> get_xs_map(boost::property_tree::ptree root);
};


#endif //READERS_H
