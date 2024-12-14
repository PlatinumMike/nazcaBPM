//
// Created by mike on 11/11/24.
//

#include "Port.h"

#include <utility>

Port::Port(std::string name, std::string placement, const double x0, std::vector<double> const &ygrid,
           std::vector<double> const &zgrid): RectangularGrid(ygrid, zgrid), name(std::move(name)), x0(x0), placement(std::move(placement)) {
}

std::string Port::get_name() const {
    return name;
}

double Port::get_x0() const {
    return x0;
}
