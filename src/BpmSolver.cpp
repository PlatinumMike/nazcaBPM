//
// Created by mike on 11/11/24.
//

#include "BpmSolver.h"
#include "FieldMonitor.h"

#include <iostream>
#include <ostream>
#include <chrono>

BpmSolver::BpmSolver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const RectangularGrid &grid,
                     const double scheme_parameter, const double k0, const double reference_index) : Solver(
    geometry, pmlx, pmly, grid, scheme_parameter, k0, reference_index, false) {
}

void BpmSolver::run(const multi_array<std::complex<double>, 2> &initial_field) const {
    auto xgrid = gridPtr->get_xgrid();
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();
    auto xgrid_m = RectangularGrid::vector_to_multi_array(xgrid);
    auto ygrid_m = RectangularGrid::vector_to_multi_array(ygrid);
    auto zgrid_m = RectangularGrid::vector_to_multi_array(zgrid);

    int numx = static_cast<int>(gridPtr->get_numx());
    int numy = static_cast<int>(gridPtr->get_numy());
    int numz = static_cast<int>(gridPtr->get_numz());

    multi_array<std::complex<double>, 2> field = initial_field;

    FieldMonitor start_field(gridPtr->get_ymin(), gridPtr->get_ymax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'x',
                             0.0,
                             numy, numz);
    start_field.populate(ygrid, zgrid, field, 0);
    start_field.save_data("field_yz_start.h5");

    //slice at z=0
    FieldMonitor field_slice_xy(gridPtr->get_xmin(), gridPtr->get_xmax(), gridPtr->get_ymin(), gridPtr->get_ymax(), 'z',
                                0.0,
                                numx, numy);
    field_slice_xy.populate(ygrid, zgrid, field, 0);

    //slice at y=0
    FieldMonitor field_slice_xz(gridPtr->get_xmin(), gridPtr->get_xmax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'y',
                                0.0,
                                numx, numz);
    field_slice_xz.populate(ygrid, zgrid, field, 0);

    const int index_1percent = numx / 100;
    const int index_10percent = numx / 10;
    const int index_50percent = numx / 2;

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int x_step = 1; x_step < numx; x_step++) {
        const double current_x = xgrid[x_step - 1];
        field = do_step_cn(field, current_x, gridPtr->get_dx());
        field_slice_xy.populate(ygrid, zgrid, field, x_step);
        field_slice_xz.populate(ygrid, zgrid, field, x_step);
        //printing rough indication of simulation progress
        if (x_step == index_1percent) {
            std::cout << "1% reached" << std::endl;
            auto end = std::chrono::steady_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            std::cout << "Elapsed time = " << delta << " (s), expected total run time = " << delta * 100 << " (s)" <<
                    std::endl;
        }
        if (x_step == index_10percent) {
            std::cout << "10% reached" << std::endl;
            auto end = std::chrono::steady_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            std::cout << "Elapsed time = " << delta << " (s), expected total run time = " << delta * 10 << " (s)" <<
                    std::endl;
        }
        if (x_step == index_50percent) {
            std::cout << "50% reached" << std::endl;
            auto end = std::chrono::steady_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            std::cout << "Elapsed time = " << delta << " (s), expected total run time = " << delta * 2 << " (s)" <<
                    std::endl;
        }
    }
    std::cout << "Engine run completed" << std::endl;
    auto end = std::chrono::steady_clock::now();
    auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
    std::cout << "Elapsed time = " << delta << " (s)" << std::endl;

    FieldMonitor end_field(gridPtr->get_ymin(), gridPtr->get_ymax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'x', 0.0,
                           numy, numz);
    end_field.populate(ygrid, zgrid, field, 0);
    end_field.save_data("field_yz_end.h5");

    field_slice_xy.save_data("field_xy.h5");
    field_slice_xz.save_data("field_xz.h5");
}
