//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "OperatorSuite.h"
#include "GridInterpolator.h"

#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

using boost::multi_array;
using boost::extents;

Solver::Solver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
               const RectangularGrid &grid,
               const double scheme_parameter, const double k0,
               const double reference_index): k0(k0),
                                              reference_index(reference_index),
                                              scheme_parameter(scheme_parameter) {
    geometryPtr = &geometry;
    pmlyPtr = &pmly;
    pmlzPtr = &pmlz;
    gridPtr = &grid;

    const int num1 = static_cast<int>(grid.get_numy());
    const int num2 = static_cast<int>(grid.get_numz());
    internal_field.resize(extents[num1][num2]);
}

// todo: error prone to skip over boundary points. off by one in the array indexing can easily occur.
// suggested solution: include the boundary points in the matrix, but overrule that row, and rhs entry to force the solution to zero on the edges.
// todo: also, a lot of code duplication here in do_step_cn, and get_rhs. Generalize, then clean up.
multi_array<std::complex<double>, 2> Solver::do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                        const double x, const double dx,
                                                        const std::complex<double> propagation_factor) const {
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();
    auto ygrid_bulk = std::vector(ygrid.begin() + 1, ygrid.end() - 1);
    auto zgrid_bulk = std::vector(zgrid.begin() + 1, zgrid.end() - 1);
    int numy = static_cast<int>(ygrid.size());
    int numz = static_cast<int>(zgrid.size());


    //cache index values
    multi_array<double, 2> index_slice_start = geometryPtr->get_index_plane(x, ygrid, zgrid);
    multi_array<double, 2> index_slice_mid = geometryPtr->get_index_plane(x + 0.5 * dx, ygrid, zgrid);
    multi_array<double, 2> index_slice_end = geometryPtr->get_index_plane(x + dx, ygrid, zgrid);

    // prefactor is p=(alpha-1)*dx*i/(2k0*n0).
    const std::complex<double> preFactorRHS = (scheme_parameter - 1.0) * propagation_factor;

    //get RHS vector, store as multi-array.
    auto rhs = OperatorSuite::get_rhs(field, ygrid, zgrid, index_slice_start, reference_index, k0, preFactorRHS,
                                      pmlyPtr,
                                      pmlzPtr);


    const auto preFactorLHS = scheme_parameter * propagation_factor;

    //solver for half step field
    multi_array<std::complex<double>, 2> half_step(extents[numy][numz]);
    std::fill_n(half_step.data(), half_step.num_elements(), 0.0);
    //invert (1+...Gy) operator, so only coupling between neighbours in y direction, this means we have numz independent problems
    for (int idz = 1; idz < numz - 1; idz++) {
        std::vector<double> position_mid(ygrid_bulk.size());
        std::vector<double> position_back(ygrid_bulk.size());
        std::vector<double> position_forward(ygrid_bulk.size());
        std::vector<double> index_mid(ygrid_bulk.size());
        std::vector<double> index_back(ygrid_bulk.size());
        std::vector<double> index_forward(ygrid_bulk.size());
        std::vector<std::complex<double> > rhs_slice(ygrid_bulk.size());
        //fill vectors
        for (auto array_index = 0; array_index < ygrid_bulk.size(); array_index++) {
            position_mid[array_index] = ygrid[array_index + 1];
            position_back[array_index] = ygrid[array_index];
            position_forward[array_index] = ygrid[array_index + 2];
            index_mid[array_index] = index_slice_mid[array_index + 1][idz];
            index_back[array_index] = index_slice_mid[array_index][idz];
            index_forward[array_index] = index_slice_mid[array_index + 2][idz];

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
        std::vector<double> position_mid(zgrid_bulk.size());
        std::vector<double> position_back(zgrid_bulk.size());
        std::vector<double> position_forward(zgrid_bulk.size());
        std::vector<double> index_mid(zgrid_bulk.size());
        std::vector<double> index_back(zgrid_bulk.size());
        std::vector<double> index_forward(zgrid_bulk.size());
        std::vector<std::complex<double> > rhs_slice(zgrid_bulk.size());
        //fill vectors
        for (auto array_index = 0; array_index < zgrid_bulk.size(); array_index++) {
            position_mid[array_index] = zgrid[array_index + 1];
            position_back[array_index] = zgrid[array_index];
            position_forward[array_index] = zgrid[array_index + 2];
            index_mid[array_index] = index_slice_end[idy][array_index + 1];
            index_back[array_index] = index_slice_end[idy][array_index];
            index_forward[array_index] = index_slice_end[idy][array_index + 2];
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

multi_array<std::complex<double>, 2> Solver::interpolate_field(const std::vector<double> &ygrid_new,
                                                               const std::vector<double> &zgrid_new) const {
    const GridInterpolator grid_interpolator(gridPtr->get_ygrid(), gridPtr->get_zgrid(), internal_field,
                                             std::complex{0.0, 0.0});
    const int num1 = static_cast<int>(ygrid_new.size());
    const int num2 = static_cast<int>(zgrid_new.size());
    multi_array<std::complex<double>, 2> new_field(extents[num1][num2]);

    for (auto i = 0; i < ygrid_new.size(); i++) {
        for (auto j = 0; j < zgrid_new.size(); j++) {
            new_field[i][j] = grid_interpolator.get_value(ygrid_new[i], zgrid_new[j]);
        }
    }

    return new_field;
}
