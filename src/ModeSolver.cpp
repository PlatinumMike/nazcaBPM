//
// Created by mike on 11/11/24.
//

#include "ModeSolver.h"

#include <iostream>
#include <format>

#include "GridInterpolator.h"

ModeSolver::ModeSolver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source,
                       const RectangularGrid &grid, const double scheme_parameter, const double k0,
                       const double reference_index, const Port &port): Solver(
                                                                            geometry, pmlx, pmly, grid,
                                                                            scheme_parameter, k0, reference_index,
                                                                            true),
                                                                        port(port), sourcePtr(&source) {
    const int num1 = static_cast<int>(grid.get_numy());
    const int num2 = static_cast<int>(grid.get_numz());
    internal_field.resize(boost::extents[num1][num2]);

    beta = 0.0;
}

void ModeSolver::run(const int max_iterations) {
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();

    std::cout << std::format("Mode solving for port {} now...", port.get_name());

    // slight offset so it excites also odd modes, this can be used to find the higher order modes
    constexpr double eps_y = 0.5;
    constexpr double eps_z = 0.3;
    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    internal_field = sourcePtr->get_initial_profile(
        ygrid, zgrid, port.get_y0() + eps_y, port.get_z0() + eps_z, 1.0, 1.0);

    for (int x_step = 0; x_step < max_iterations; x_step++) {
        internal_field = do_step_cn(internal_field, port.get_x0(), gridPtr->get_dx());
        //todo: break from loop if beta is in steady state
        // beta = ...
    }
}

double ModeSolver::get_beta() const {
    return beta;
}

multi_array<std::complex<double>, 2> ModeSolver::interpolate_field(const std::vector<double> &ygrid,
                                                                   const std::vector<double> &zgrid) const {
    const GridInterpolator grid_interpolator(gridPtr->get_ygrid(), gridPtr->get_zgrid(), internal_field,
                                             std::complex{0.0, 0.0});
    multi_array<std::complex<double>, 2> new_field;

    for (auto i = 0; i < ygrid.size(); i++) {
        for (auto j = 0; j < zgrid.size(); j++) {
            new_field[i][j] = grid_interpolator.get_value(ygrid[i], zgrid[j]);
        }
    }

    return new_field;
}
