//
// Created by mike on 10/11/24.
//

#ifndef MODEHANDLER_H
#define MODEHANDLER_H

#include <boost/multi_array.hpp>
using boost::multi_array;

class ModeHandler {
public:
    explicit ModeHandler(const std::string &source_type);

    /**
     * Generates a 2D Gaussian profile
     * @param xgrid grid of x positions
     * @param ygrid grid of y positions
     * @param x0 center of Gaussian in x
     * @param y0 center of Gaussian in y
     * @param std_x standard deviation in x
     * @param std_y standard deviation in y
     * @return
     */
    multi_array<std::complex<double>, 2> get_initial_profile(const std::vector<double> &xgrid,
                                                             const std::vector<double> &ygrid, double x0,
                                                             double y0, double std_x,
                                                             double std_y) const;

private:
    std::string type;
    //todo: add option to read a mode profile from a file, then interpolate on the BPM mesh.
    //todo: add mode overlap integrals
};


#endif //MODEHANDLER_H
