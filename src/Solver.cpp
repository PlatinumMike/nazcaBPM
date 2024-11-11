//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "OperatorSuite.h"

#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;

Solver::Solver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
               const RectangularGrid &grid,
               const double scheme_parameter, const double k0,
               const double reference_index, const bool mode_solve): scheme_parameter(scheme_parameter),
                                                                     k0(k0), reference_index(reference_index),
                                                                     mode_solve(mode_solve) {
    geometryPtr = &geometry;
    pmlyPtr = &pmly;
    pmlzPtr = &pmlz;
    gridPtr = &grid;
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

    std::complex<double> base = dx / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};
    if (mode_solve) {
        //For the Imaginary distance method we propagate along ix instead of x.
        base *= std::complex<double>{0.0, 1.0};
    }

    // get index in slice
    multi_array<double, 2> index_slice(extents[numy][numz]);
    for (int idy = 0; idy < numy; idy++) {
        for (int idz = 0; idz < numz; idz++) {
            index_slice[idy][idz] = geometryPtr->get_index(x, ygrid[idy], zgrid[idz]);
        }
    }
    // prefactor is p=(alpha-1)*dx*i/(2k0*n0).
    const std::complex<double> preFactorRHS = (scheme_parameter - 1.0) * base;

    //get RHS vector, store as multi-array.
    auto rhs = OperatorSuite::get_rhs(field, ygrid, zgrid, index_slice, reference_index, k0, preFactorRHS, pmlyPtr,
                                      pmlzPtr);


    const auto preFactorLHS = scheme_parameter * base;

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
