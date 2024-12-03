//
// Created by mike on 11/11/24.
//

#include "BpmSolver.h"
#include "FieldMonitor.h"
#include "ProgressBar.h"

BpmSolver::BpmSolver(const Geometry &geometry, const PML &pmly, const PML &pmlz, const RectangularGrid3D &grid,
                     const double scheme_parameter, const double k0, const double reference_index,
                     const std::filesystem::path &absolute_path_output) : Solver(
                                                                              geometry, pmly, pmlz, grid,
                                                                              scheme_parameter, k0,
                                                                              reference_index), grid3d(grid),
                                                                          absolute_path_output(absolute_path_output) {
}

void BpmSolver::run(const multi_array<std::complex<double>, 2> &initial_field) {
    auto xgrid = grid3d.get_xgrid();
    auto ygrid = grid3d.get_ygrid();
    auto zgrid = grid3d.get_zgrid();

    int numx = static_cast<int>(grid3d.get_numx());
    int numy = static_cast<int>(grid3d.get_numy());
    int numz = static_cast<int>(grid3d.get_numz());

    std::filesystem::path output_path = absolute_path_output;

    internal_field = initial_field;

    FieldMonitor start_field(grid3d.get_ymin(), grid3d.get_ymax(), grid3d.get_zmin(), grid3d.get_zmax(), 'x',
                             0.0,
                             numy, numz);
    start_field.populate(ygrid, zgrid, internal_field, 0);
    output_path.append("field_yz_start.h5");
    start_field.save_data(output_path.string());

    //slice at z=0
    FieldMonitor field_slice_xy(grid3d.get_xmin(), grid3d.get_xmax(), grid3d.get_ymin(), grid3d.get_ymax(), 'z',
                                0.0,
                                numx, numy);
    field_slice_xy.populate(ygrid, zgrid, internal_field, 0);

    //slice at y=0
    FieldMonitor field_slice_xz(grid3d.get_xmin(), grid3d.get_xmax(), grid3d.get_zmin(), grid3d.get_zmax(), 'y',
                                0.0,
                                numx, numz);
    field_slice_xz.populate(ygrid, zgrid, internal_field, 0);


    const std::complex<double> propagation_factor = grid3d.get_dx() / (2.0 * k0 * reference_index)
                                                    * std::complex<double>{0.0, 1.0};

    ProgressBar progress_bar(numx);
    for (int x_step = 1; x_step < numx; x_step++) {
        const double current_x = xgrid[x_step - 1];
        internal_field = do_step_cn(internal_field, current_x, grid3d.get_dx(), propagation_factor);
        field_slice_xy.populate(ygrid, zgrid, internal_field, x_step);
        field_slice_xz.populate(ygrid, zgrid, internal_field, x_step);
        progress_bar.update(x_step);
    }
    progress_bar.finalize();

    FieldMonitor end_field(grid3d.get_ymin(), grid3d.get_ymax(), grid3d.get_zmin(), grid3d.get_zmax(), 'x', 0.0,
                           numy, numz);
    end_field.populate(ygrid, zgrid, internal_field, 0);
    output_path.replace_filename("field_yz_end.h5");
    end_field.save_data(output_path.string());

    output_path.replace_filename("field_xy.h5");
    field_slice_xy.save_data(output_path.string());
    output_path.replace_filename("field_xz.h5");
    field_slice_xz.save_data(output_path.string());
}
