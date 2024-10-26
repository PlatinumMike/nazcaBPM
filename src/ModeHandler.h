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
     * @param ygrid grid of y positions
     * @param zgrid grid of z positions
     * @param y0 center of Gaussian in y
     * @param z0 center of Gaussian in z
     * @param std_y standard deviation in y
     * @param std_z standard deviation in z
     * @return
     */
    multi_array<std::complex<double>, 2> get_initial_profile(const std::vector<double> &ygrid,
                                                             const std::vector<double> &zgrid, double y0,
                                                             double z0, double std_y,
                                                             double std_z) const;

private:
    std::string type;
    //todo: add option to read a mode profile from a file, then interpolate on the BPM mesh.
    //todo: add mode overlap integrals
};


#endif //MODEHANDLER_H
