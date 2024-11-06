//
// Created by mike on 11/5/24.
//

#ifndef FIELDMONITOR_H
#define FIELDMONITOR_H
#include <vector>
#include <boost/multi_array.hpp>

//TODO: this differs from the index monitor in that we cannot evaluate it at any point.
// it has to be on the simulation xgrid. In principle we could interpolate, but that makes things a lot more complicated.
// to be added later on. Well, in the populate function you can pass in the current x value, so that will be same for the grid an monitor.
// in the other two directions we can easily interpolate.
// Currently no check is done to see if the xgrid is the same as grid1 if orientation='x'. Better would be pass the whole grid as input.
// then we are sure that the xgrids are identical.

// TODO: for ports do this:
// Also, it needs to have way to normalize the field, such that <u,u> = 1. void normalize_field();
// And it should be able to compute modes.
// Actually, this does not need to be part of the field monitor, but it should be for a port.

// TODO: lots of duplicate code with the index monitor, use inheritance?


class FieldMonitor {
public:
    FieldMonitor(double coordinate1min, double coordinate1max, double coordinate2min, double coordinate2max,
                 char orientation, double slice_position, int resolution1, int resolution2);

    /**
     * You do not want to store then entire 3D grid in memory, so you have to construct the field of the monitor
     * from yz slices of the field that you pass in.
     *
     * if the orientation is 'x', an extra argument idx is needed to determine if the array is fully populated yet.
     */
    void populate(const std::vector<double> &ygrid, const std::vector<double> &zgrid,
                  const boost::multi_array<std::complex<double>, 2> &field_slice, int idx);

    void save_data(const std::string &filename) const;

private:
    const double coordinate1min;
    const double coordinate1max;
    const double coordinate1span;
    const double coordinate2min;
    const double coordinate2max;
    const double coordinate2span;
    const double slice_position;
    const char orientation;
    bool populated = false;

    std::vector<double> grid1;
    std::vector<double> grid2;
    std::string label1 = "grid1";
    std::string label2 = "grid2";

    boost::multi_array<std::complex<double>, 2> field;
};


#endif //FIELDMONITOR_H
