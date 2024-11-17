//
// Created by mike on 11/11/24.
//

#include "Port.h"

Port::Port(const std::string &name, const std::string &placement, const double x0, std::vector<double> const &ygrid,
           std::vector<double> const &zgrid): RectangularGrid(ygrid, zgrid), name(name), x0(x0), placement(placement) {
}

std::string Port::get_name() const {
    return name;
}

double Port::get_x0() const {
    return x0;
}
