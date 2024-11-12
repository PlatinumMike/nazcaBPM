//
// Created by mike on 11/11/24.
//

#include "Port.h"

Port::Port(const std::string &name, const std::string &placement, const double x0, const double y0, const double z0,
           const double yspan,
           const double zspan): name(name), x0(x0),
                                y0(y0), z0(z0), yspan(yspan), zspan(zspan), placement(placement) {
}

std::string Port::get_name() const {
    return name;
}

double Port::get_x0() const {
    return x0;
}

double Port::get_y0() const {
    return y0;
}

double Port::get_yspan() const {
    return yspan;
}

double Port::get_ymin() const {
    return y0 - 0.5 * yspan;
}

double Port::get_ymax() const {
    return y0 + 0.5 * yspan;
}

double Port::get_zmin() const {
    return z0 - 0.5 * zspan;
}

double Port::get_zmax() const {
    return z0 + 0.5 * zspan;
}

double Port::get_z0() const {
    return z0;
}

double Port::get_zspan() const {
    return zspan;
}
