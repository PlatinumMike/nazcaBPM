//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "OperatorSuite.h"
#include "FieldMonitor.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <chrono>

#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;

Solver::Solver(const Geometry &geometry, const PML &pmly, const PML &pmlz, const ModeHandler &source,
               const RectangularGrid &grid,
               const double scheme_parameter, const double k0,
               const double reference_index): scheme_parameter(scheme_parameter),
                                              k0(k0), reference_index(reference_index) {
    geometryPtr = &geometry;
    pmlyPtr = &pmly;
    pmlzPtr = &pmlz;
    sourcePtr = &source;
    gridPtr = &grid;
}

void Solver::run() {
    auto xgrid = gridPtr->get_xgrid();
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();
    auto xgrid_m = RectangularGrid::vector_to_multi_array(xgrid);
    auto ygrid_m = RectangularGrid::vector_to_multi_array(ygrid);
    auto zgrid_m = RectangularGrid::vector_to_multi_array(zgrid);

    int numx = static_cast<int>(gridPtr->get_numx());
    int numy = static_cast<int>(gridPtr->get_numy());
    int numz = static_cast<int>(gridPtr->get_numz());

    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    multi_array<std::complex<double>, 2> field = sourcePtr->get_initial_profile(ygrid, zgrid, 0.0, 0.0, 1.0, 1.0);

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

// todo: error prone to skip over boundary points. off by one in the array indexing can easily occur.
// suggested solution: include the boundary points in the matrix, but overrule that row, and rhs entry to force the solution to zero on the edges.
// todo: also, a lot of code duplication here in do_step_cn, and get_rhs. Generalize, then clean up.
multi_array<std::complex<double>, 2> Solver::do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                        const double x, const double dx) const {
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();
    int numy = static_cast<int>(ygrid.size());
    int numz = static_cast<int>(zgrid.size());

    // get index in slice
    multi_array<double, 2> index_slice(extents[numy][numz]);
    for (int idy = 0; idy < numy; idy++) {
        for (int idz = 0; idz < numz; idz++) {
            index_slice[idy][idz] = geometryPtr->get_index(x, ygrid[idy], zgrid[idz]);
        }
    }
    // prefactor is p=(alpha-1)*dx*i/(2k0*n0).
    const std::complex<double> preFactorRHS =
            (scheme_parameter - 1.0) * dx / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //get RHS vector, store as multi-array.
    auto rhs = OperatorSuite::get_rhs(field, ygrid, zgrid, index_slice, reference_index, k0, preFactorRHS, pmlyPtr,
                                      pmlzPtr);


    const auto preFactorLHS = scheme_parameter * dx / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //solver for half step field
    multi_array<std::complex<double>, 2> half_step(extents[numy][numz]);
    std::fill_n(half_step.data(), half_step.num_elements(), 0.0);
    //invert (1+...Gy) operator, so only coupling between neighbours in y direction, this means we have numz independent problems
    for (int idz = 1; idz < numz - 1; idz++) {
        auto size = numy - 2;
        std::vector<double> position_mid(size);
        std::vector<double> position_back(size);
        std::vector<double> position_forward(size);
        std::vector<double> index_mid(size);
        std::vector<double> index_back(size);
        std::vector<double> index_forward(size);
        std::vector<std::complex<double> > rhs_slice(size);
        //fill vectors
        for (auto array_index = 0; array_index < size; array_index++) {
            position_mid[array_index] = ygrid[array_index + 1];
            position_back[array_index] = ygrid[array_index];
            position_forward[array_index] = ygrid[array_index + 2];
            index_mid[array_index] = geometryPtr->get_index(x + 0.5 * dx, ygrid[array_index + 1], zgrid[idz]);
            index_back[array_index] = geometryPtr->get_index(x + 0.5 * dx, ygrid[array_index], zgrid[idz]);
            index_forward[array_index] = geometryPtr->get_index(x + 0.5 * dx, ygrid[array_index + 2], zgrid[idz]);

            rhs_slice[array_index] = rhs[array_index + 1][idz];
        }
        auto solution = OperatorSuite::solve_system(position_mid, position_back, position_forward, index_mid,
                                                    index_back, index_forward, reference_index, k0, preFactorLHS,
                                                    rhs_slice, pmlyPtr);


        //store into solution vector
        for (int idy = 1; idy < numy - 1; idy++) {
            half_step[idy][idz] = solution[idy - 1];
        }
    }


    // solve for next field value.
    multi_array<std::complex<double>, 2> full_step(extents[numy][numz]);
    std::fill_n(full_step.data(), full_step.num_elements(), 0.0);
    //invert (1+...Gz) operator, so only coupling between neighbours in z direction, this means we have numy independent problems
    for (int idy = 1; idy < numy - 1; idy++) {
        auto size = numz - 2;
        std::vector<double> position_mid(size);
        std::vector<double> position_back(size);
        std::vector<double> position_forward(size);
        std::vector<double> index_mid(size);
        std::vector<double> index_back(size);
        std::vector<double> index_forward(size);
        std::vector<std::complex<double> > rhs_slice(size);
        //fill vectors
        for (auto array_index = 0; array_index < size; array_index++) {
            position_mid[array_index] = zgrid[array_index + 1];
            position_back[array_index] = zgrid[array_index];
            position_forward[array_index] = zgrid[array_index + 2];
            index_mid[array_index] = geometryPtr->get_index(x + dx, ygrid[idy], zgrid[array_index + 1]);
            index_back[array_index] = geometryPtr->get_index(x + dx, ygrid[idy], zgrid[array_index]);
            index_forward[array_index] = geometryPtr->get_index(x + dx, ygrid[idy], zgrid[array_index + 2]);
            //using the half_step now to fill the rhs vector.
            rhs_slice[array_index] = half_step[idy][array_index + 1];
        }

        auto solution = OperatorSuite::solve_system(position_mid, position_back, position_forward, index_mid,
                                                    index_back, index_forward, reference_index, k0, preFactorLHS,
                                                    rhs_slice, pmlzPtr);


        //store into solution vector
        for (int idz = 1; idz < numz - 1; idz++) {
            full_step[idy][idz] = solution[idz - 1];
        }
    }
    return full_step;
}


multi_array<double, 2> Solver::get_intensity(const multi_array<std::complex<double>, 2> &field) const {
    const int numy = static_cast<int>(field.shape()[0]);
    const int numz = static_cast<int>(field.shape()[1]);

    multi_array<double, 2> intensity(extents[numy][numz]);
    double max = 0.0;
    for (auto idy = 0; idy < field.shape()[0]; idy++) {
        for (auto idz = 0; idz < field.shape()[1]; idz++) {
            const double real = field[idy][idz].real();
            const double imag = field[idy][idz].imag();
            const double value = real * real + imag * imag;
            intensity[idy][idz] = value;
            max = std::max(max, value);
        }
    }
    for (auto idy = 0; idy < field.shape()[0]; idy++) {
        for (auto idz = 0; idz < field.shape()[1]; idz++) {
            intensity[idy][idz] /= max;
        }
    }
    return intensity;
}
