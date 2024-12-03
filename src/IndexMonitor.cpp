//
// Created by mike on 11/5/24.
//

#include "IndexMonitor.h"
#include "AuxiliaryFunctions.h"
#include "RectangularGrid.h"
#include "IO/hdf_writer.h"

#include <cassert>


IndexMonitor::IndexMonitor(const double coordinate1min, const double coordinate1max, const double coordinate2min,
                           const double coordinate2max, const char orientation, const double slice_position,
                           const int resolution1,
                           const int resolution2): coordinate1min(coordinate1min),
                                                   coordinate1max(coordinate1max),
                                                   coordinate1span(coordinate1max - coordinate1min),
                                                   coordinate2min(coordinate2min),
                                                   coordinate2max(coordinate2max),
                                                   coordinate2span(coordinate2max - coordinate2min),
                                                   slice_position(slice_position),
                                                   orientation(orientation) {
    // monitor must have positive sizes
    assert(coordinate1span > 0.0);
    assert(coordinate2span > 0.0);
    assert(resolution1>0);
    assert(resolution2>0);
    assert(
        ("Monitor orientation not recognized, use x or y or z", orientation=='x' || orientation=='y'|| orientation=='z'
        ));

    grid1 = AuxiliaryFunctions::linspace(coordinate1min, coordinate1max, resolution1);
    grid2 = AuxiliaryFunctions::linspace(coordinate2min, coordinate2max, resolution2);

    if (orientation == 'x') {
        label1 = "ygrid";
        label2 = "zgrid";
    } else if (orientation == 'y') {
        label1 = "xgrid";
        label2 = "zgrid";
    } else {
        label1 = "xgrid";
        label2 = "ygrid";
    }

    const int num1 = static_cast<int>(grid1.size());
    const int num2 = static_cast<int>(grid2.size());
    index_dataset.resize(boost::extents[num1][num2]);
}

void IndexMonitor::populate(const Geometry &geometryPtr) {
    if (populated) {
        printf("Index monitor already populated! Skipping\n");
        return;
    }

    for (auto index1 = 0; index1 < index_dataset.shape()[0]; index1++) {
        for (auto index2 = 0; index2 < index_dataset.shape()[1]; index2++) {
            if (orientation == 'x') {
                index_dataset[index1][index2] = geometryPtr.get_index(slice_position, grid1[index1],
                                                                      grid2[index2]);
            } else if (orientation == 'y') {
                index_dataset[index1][index2] = geometryPtr.get_index(grid1[index1], slice_position,
                                                                      grid2[index2]);
            } else {
                index_dataset[index1][index2] = geometryPtr.get_index(grid1[index1],
                                                                      grid2[index2],
                                                                      slice_position);
            }
        }
    }
    populated = true;
}

double IndexMonitor::get_max_index() {
    double max = index_dataset[0][0];
    for (auto id1 = 0; id1 < index_dataset.shape()[0]; id1++) {
        for (auto id2 = 0; id2 < index_dataset.shape()[1]; id2++) {
            max = std::max(max, index_dataset[id1][id2]);
        }
    }
    return max;
}

double IndexMonitor::get_min_index() {
    double min = index_dataset[0][0];
    for (auto id1 = 0; id1 < index_dataset.shape()[0]; id1++) {
        for (auto id2 = 0; id2 < index_dataset.shape()[1]; id2++) {
            min = std::min(min, index_dataset[id1][id2]);
        }
    }
    return min;
}

void IndexMonitor::save_data(const std::string &filename) const {
    if (!(populated)) {
        printf("Index monitor not yet populated! Cannot save. Skipping\n");
        return;
    }
    const auto grid1_m = RectangularGrid::vector_to_multi_array(grid1);
    const auto grid2_m = RectangularGrid::vector_to_multi_array(grid2);
    H5::H5File file(filename, H5F_ACC_TRUNC);
    write_hdf5(file, "refractive_index", index_dataset);
    write_hdf5(file, label1, grid1_m);
    write_hdf5(file, label2, grid2_m);
    file.close();
}
