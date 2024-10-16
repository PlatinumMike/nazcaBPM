//
// Created by mike on 10/11/24.
//

#include "ModeHandler.h"

#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;


ModeHandler::ModeHandler(const std::string &source_type): type(source_type) {
}

// just some mock-up profile.
multi_array<std::complex<double>, 2> ModeHandler::get_initial_profile(const std::vector<double> &xgrid,
                                                                      const std::vector<double> &ygrid, const double x0,
                                                                      const double y0, const double std_x,
                                                                      const double std_y) const {
    //todo: in the future this would be read from a file, which is generated by a mode solver, e.g. WGMS3D.
    //still it is convenient to keep to keep this mockup profile. Maybe handy as initial profile for mode finding (imaginary distance BPM).

    int numx = static_cast<int>(xgrid.size());
    int numy = static_cast<int>(ygrid.size());
    multi_array<std::complex<double>, 2> initial_profile(extents[numx][numy]);
    std::complex<double> value = {0.0, 0.0};
    double x = 0.0;
    double y = 0.0;
    for (int idx = 0; idx < numx; idx++) {
        for (int idy = 0; idy < numy; idy++) {
            if (idx == 0 || idx == numx - 1 || idy == 0 || idy == numy - 1) {
                //set boundary values to zero. Simple Dirichlet BCs. Reflections are negligible since a PML will be placed in front of the metal wall anyway.
                value = {0.0, 0.0};
            } else {
                x = xgrid[idx];
                y = ygrid[idy];
                double deltax = (x - x0) / std_x;
                double deltay = (y - y0) / std_y;
                value = {std::exp(-0.5 * deltax * deltax - 0.5 * deltay * deltay), 0.0};
            }
            initial_profile[idx][idy] = value;
        }
    }
    return initial_profile;
}