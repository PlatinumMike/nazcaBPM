//
// Created by mike on 11/5/24.
//

#ifndef INDEXMONITOR_H
#define INDEXMONITOR_H

#include <string>
#include <vector>
#include <boost/multi_array.hpp>


#include "Geometry.h"

class IndexMonitor {
public:
    /**
     * Define an index monitor. This is a flat region in which the index is computed for later analysis.
     * It creates its own grid, so in principle that can be any resolution. However, it is recommended to use the same
     * resolution as in the rest of the simulation so you get a granularity that is representative of what the actual
     * simulation uses under the hood.
     *
     * The coordinates can be x,y or z. The orientation must be along one of the unit vectors.
     * The slice_position indicates where the monitor is placed along the axis that is chosen by the "orientation".
     */
    IndexMonitor(double coordinate1min, double coordinate1max, double coordinate2min, double coordinate2max,
                 char orientation, double slice_position, int resolution1, int resolution2);

    void populate(const Geometry &geometryPtr);

    double get_max_index();

    double get_min_index();

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

    boost::multi_array<double, 2> index_dataset;
};


#endif //INDEXMONITOR_H
