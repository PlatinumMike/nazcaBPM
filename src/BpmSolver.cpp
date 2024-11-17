//
// Created by mike on 11/11/24.
//

#include "BpmSolver.h"
#include "FieldMonitor.h"
#include "ProgressBar.h"

BpmSolver::BpmSolver(const Geometry &geometry, const PML &pmly, const PML &pmlz, const RectangularGrid &grid,
                     const double scheme_parameter, const double k0, const double reference_index) : Solver(
    geometry, pmly, pmlz, grid, scheme_parameter, k0, reference_index) {
}

void BpmSolver::run(const multi_array<std::complex<double>, 2> &initial_field) {
    auto xgrid = gridPtr->get_xgrid();
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();
    auto xgrid_m = RectangularGrid::vector_to_multi_array(xgrid);
    auto ygrid_m = RectangularGrid::vector_to_multi_array(ygrid);
    auto zgrid_m = RectangularGrid::vector_to_multi_array(zgrid);

    int numx = static_cast<int>(gridPtr->get_numx());
    int numy = static_cast<int>(gridPtr->get_numy());
    int numz = static_cast<int>(gridPtr->get_numz());

    internal_field = initial_field;

    FieldMonitor start_field(gridPtr->get_ymin(), gridPtr->get_ymax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'x',
                             0.0,
                             numy, numz);
    start_field.populate(ygrid, zgrid, internal_field, 0);
    start_field.save_data("field_yz_start.h5");

    //slice at z=0
    FieldMonitor field_slice_xy(gridPtr->get_xmin(), gridPtr->get_xmax(), gridPtr->get_ymin(), gridPtr->get_ymax(), 'z',
                                0.0,
                                numx, numy);
    field_slice_xy.populate(ygrid, zgrid, internal_field, 0);

    //slice at y=0
    FieldMonitor field_slice_xz(gridPtr->get_xmin(), gridPtr->get_xmax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'y',
                                0.0,
                                numx, numz);
    field_slice_xz.populate(ygrid, zgrid, internal_field, 0);


    const std::complex<double> propagation_factor = gridPtr->get_dx() / (2.0 * k0 * reference_index)
                                                    * std::complex<double>{0.0, 1.0};

    ProgressBar progress_bar(numx);
    for (int x_step = 1; x_step < numx; x_step++) {
        const double current_x = xgrid[x_step - 1];
        internal_field = do_step_cn(internal_field, current_x, gridPtr->get_dx(), propagation_factor);
        field_slice_xy.populate(ygrid, zgrid, internal_field, x_step);
        field_slice_xz.populate(ygrid, zgrid, internal_field, x_step);
        progress_bar.update(x_step);
    }
    progress_bar.finalize();

    FieldMonitor end_field(gridPtr->get_ymin(), gridPtr->get_ymax(), gridPtr->get_zmin(), gridPtr->get_zmax(), 'x', 0.0,
                           numy, numz);
    end_field.populate(ygrid, zgrid, internal_field, 0);
    end_field.save_data("field_yz_end.h5");

    field_slice_xy.save_data("field_xy.h5");
    field_slice_xz.save_data("field_xz.h5");
}
