//
// Created by mike on 11/5/24.
//

#include "FieldMonitor.h"
#include "AuxiliaryFunctions.h"
#include "GridInterpolator.h"
#include "RectangularGrid.h"
#include "IO/hdf_writer.h"

FieldMonitor::FieldMonitor(const double coordinate1min, const double coordinate1max, const double coordinate2min,
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
    field.resize(boost::extents[num1][num2]);
}

void FieldMonitor::populate(const std::vector<double> &ygrid, const std::vector<double> &zgrid,
                            const multi_array<std::complex<double>, 2> &field_slice, const int idx) {
    if (populated) {
        printf("Field monitor already populated! Skipping\n");
        return;
    }

    // note, the monitor can be oriented any direction, but the field_slice is defined on the y,z grid.
    // so when interpolating it you must pass in an y,z value.
    const GridInterpolator grid_interpolator(ygrid, zgrid, field_slice, std::complex{0.0, 0.0});
    if (orientation == 'x') {
        //yz monitor, so grid1 = ygrid, grid2 = zgrid
        //slice_position is irrelevant because we pass in a whole slice. So the x position is whatever we are currently at.
        for (auto i = 0; i < field.shape()[0]; i++) {
            for (auto j = 0; j < field.shape()[1]; j++) {
                field[i][j] = grid_interpolator.get_value(grid1[i], grid2[j]);
            }
        }
        populated = true;
    } else if (orientation == 'y') {
        //xz monitor, so grid1 = xgrid, grid2 = zgrid
        //only loop over second coordinate because we just have one x value
        for (auto j = 0; j < field.shape()[1]; j++) {
            field[idx][j] = grid_interpolator.get_value(slice_position, grid2[j]);
        }
        populated = idx == static_cast<int>(grid1.size() - 1);
    } else {
        //xy monitor, so grid1 = xgrid, grid2 = ygrid
        for (auto j = 0; j < field.shape()[1]; j++) {
            field[idx][j] = grid_interpolator.get_value(grid2[j], slice_position);
        }
        populated = idx == static_cast<int>(grid1.size() - 1);
    }
}


void FieldMonitor::save_data(const std::string &filename) const {
    if (!(populated)) {
        printf("Field monitor not yet populated! Cannot save. Skipping\n");
        return;
    }
    const auto grid1_m = RectangularGrid::vector_to_multi_array(grid1);
    const auto grid2_m = RectangularGrid::vector_to_multi_array(grid2);
    write_cmplx_hdf5(filename, field, grid1_m, grid2_m, label1, label2);
}
