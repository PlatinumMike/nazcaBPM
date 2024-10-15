//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "AuxiliaryFunctions.h"
#include "IO/Readers.h"
#include "IO/hdf_writer.h"
#include "OperatorSuite.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <chrono>

#include <hdf5/serial/H5Cpp.h>
#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;

Solver::Solver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source, double xmin,
               double xmax, double ymin, double ymax, double zmin, double zmax, int numx, int numy,
               int numz, const double scheme_parameter, double k0, double reference_index): xmin(xmin), xmax(xmax),
    ymin(ymin), ymax(ymax), zmin(zmin),
    zmax(zmax), numx(numx),
    numy(numy), numz(numz), scheme_parameter(scheme_parameter), k0(k0), reference_index(reference_index) {
    geometryPtr = &geometry;
    pmlxPtr = &pmlx;
    pmlyPtr = &pmly;
    sourcePtr = &source;
}

void Solver::run() {
    auto xgrid = AuxiliaryFunctions::linspace(xmin, xmax, numx);
    auto ygrid = AuxiliaryFunctions::linspace(ymin, ymax, numy);
    auto zgrid = AuxiliaryFunctions::linspace(zmin, zmax, numz);
    dx = xgrid[1] - xgrid[0];
    dy = ygrid[1] - ygrid[0];
    dz = zgrid[1] - zgrid[0];

    auto xgrid_m = vector_to_multi_array(xgrid);
    auto ygrid_m = vector_to_multi_array(ygrid);
    auto zgrid_m = vector_to_multi_array(zgrid);

    dump_index_slice("index_start.h5", 'z', 0.0, xgrid, ygrid);
    dump_index_slice("index_end.h5", 'z', zgrid[numz - 1], xgrid, ygrid);
    dump_index_slice("index_yz.h5", 'x', 0.0, ygrid, zgrid);
    dump_index_slice("index_xz.h5", 'y', 0.0, xgrid, zgrid);

    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    multi_array<std::complex<double>, 2> field = sourcePtr->get_initial_profile(xgrid, ygrid, 0.0, 0.0, 1.0, 1.0);

    write_cmplx_hdf5("field_start.h5", field, xgrid_m, ygrid_m, 'z');

    multi_array<std::complex<double>, 2> field_slice_yz(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_slice_xz(extents[numx][numz]);
    record_slice(field, field_slice_yz, 0, true);
    record_slice(field, field_slice_xz, 0, false);


    const int index_1percent = numz / 100;
    const int index_10percent = numz / 10;
    const int index_50percent = numz / 2;

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int z_step = 1; z_step < numz; z_step++) {
        const double current_z = zgrid[z_step - 1];
        field = do_step_cn(field, xgrid, ygrid, current_z, dz);
        record_slice(field, field_slice_yz, z_step, true);
        record_slice(field, field_slice_xz, z_step, false);
        //printing rough indication of simulation progress
        if (z_step == index_1percent) {
            std::cout << "1% reached" << std::endl;
            auto end = std::chrono::steady_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            std::cout << "Elapsed time = " << delta << " (s), expected total run time = " << delta * 100 << " (s)" <<
                    std::endl;
        }
        if (z_step == index_10percent) {
            std::cout << "10% reached" << std::endl;
            auto end = std::chrono::steady_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            std::cout << "Elapsed time = " << delta << " (s), expected total run time = " << delta * 10 << " (s)" <<
                    std::endl;
        }
        if (z_step == index_50percent) {
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

    write_cmplx_hdf5("field_end.h5", field, xgrid_m, ygrid_m, 'z');
    write_cmplx_hdf5("field_xz.h5", field_slice_xz, xgrid_m, zgrid_m, 'y');
    write_cmplx_hdf5("field_yz.h5", field_slice_yz, ygrid_m, zgrid_m, 'x');
}

double Solver::get_refractive_index(const double x_bpm, const double y_bpm, double z_bpm) const {
    //convert to nazca coordinates
    return geometryPtr->get_index(z_bpm, -y_bpm, x_bpm);
}


void Solver::record_slice(const multi_array<std::complex<double>, 2> &buffer,
                          multi_array<std::complex<double>, 2> &storage, const int idz, const bool slice_x) const {
    if (slice_x) {
        const int index_xmid = numx / 2; //yz slice at x=0 approximately
        for (int idy = 0; idy < numy; idy++) {
            storage[idy][idz] = buffer[index_xmid][idy];
        }
    } else {
        const int index_ymid = numy / 2; //xz slice at y=0 approximately
        for (int idx = 0; idx < numx; idx++) {
            storage[idx][idz] = buffer[idx][index_ymid];
        }
    }
}

void Solver::dump_index_slice(const std::string &filename, const char direction, const double slice_position,
                              const std::vector<double> &grid_coordinate1,
                              const std::vector<double> &grid_coordinate2) const {
    assert(("Slice direction not recognized, use x or y or z", direction=='x' || direction=='y'|| direction=='z'));

    const int num1 = static_cast<int>(grid_coordinate1.size());
    const int num2 = static_cast<int>(grid_coordinate2.size());
    multi_array<double, 2> index_dataset(extents[num1][num2]);
    for (int index1 = 0; index1 < num1; index1++) {
        for (int index2 = 0; index2 < num2; index2++) {
            if (direction == 'x') {
                index_dataset[index1][index2] = get_refractive_index(slice_position, grid_coordinate1[index1],
                                                                     grid_coordinate2[index2]);
            } else if (direction == 'y') {
                index_dataset[index1][index2] = get_refractive_index(grid_coordinate1[index1], slice_position,
                                                                     grid_coordinate2[index2]);
            } else {
                index_dataset[index1][index2] = get_refractive_index(grid_coordinate1[index1], grid_coordinate2[index2],
                                                                     slice_position);
            }
        }
    }

    std::string label1 = "grid1";
    std::string label2 = "grid2";
    if (direction == 'x') {
        label1 = "ygrid";
        label2 = "zgrid";
    } else if (direction == 'y') {
        label1 = "xgrid";
        label2 = "zgrid";
    } else {
        label1 = "xgrid";
        label2 = "ygrid";
    }

    const auto grid1 = vector_to_multi_array(grid_coordinate1);
    const auto grid2 = vector_to_multi_array(grid_coordinate2);
    H5::H5File file(filename, H5F_ACC_TRUNC);
    write_hdf5(file, "refractive_index", index_dataset);
    write_hdf5(file, label1, grid1);
    write_hdf5(file, label2, grid2);
    file.close();
}

multi_array<double, 1> Solver::vector_to_multi_array(const std::vector<double> &vec) {
    const int size = static_cast<int>(vec.size());
    multi_array<double, 1> vec_m(extents[size]);
    for (int i = 0; i < size; i++) {
        vec_m[i] = vec[i];
    }
    return vec_m;
}

// todo: error prone to skip over boundary points. off by one in the array indexing can easily occur.
// suggested solution: include the boundary points in the matrix, but overrule that row, and rhs entry to force the solution to zero on the edges.
// todo: also, a lot of code duplication here in do_step_cn, and get_rhs. Generalize, then clean up.
multi_array<std::complex<double>, 2> Solver::do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                        const std::vector<double> &xgrid,
                                                        const std::vector<double> &ygrid, const double z,
                                                        const double dz) const {
    // get index in slice
    multi_array<double, 2> index_slice(extents[numx][numy]);
    for (int idx = 0; idx < numx; idx++) {
        for (int idy = 0; idy < numy; idy++) {
            index_slice[idx][idy] = get_refractive_index(xgrid[idx], ygrid[idy], z);
        }
    }
    // prefactor is p=(alpha-1)*dz*i/(2k0*n0).
    const std::complex<double> preFactorRHS =
            (scheme_parameter - 1.0) * dz / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //get RHS vector, store as multi-array.
    auto rhs = OperatorSuite::get_rhs(field, xgrid, ygrid, index_slice, reference_index, k0, preFactorRHS, pmlxPtr,
                                      pmlyPtr);


    const auto preFactorLHS = scheme_parameter * dz / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //solver for half step field
    multi_array<std::complex<double>, 2> half_step(extents[numx][numy]);
    std::fill_n(half_step.data(), half_step.num_elements(), 0.0);
    //invert (1+...Gx) operator, so only coupling between neighbours in x direction, this means we have numy independent problems
    for (int idy = 1; idy < numy - 1; idy++) {
        auto size = numx - 2;
        std::vector<double> position_mid(size);
        std::vector<double> position_back(size);
        std::vector<double> position_forward(size);
        std::vector<double> index_mid(size);
        std::vector<double> index_back(size);
        std::vector<double> index_forward(size);
        std::vector<std::complex<double> > rhs_slice(size);
        //fill vectors
        for (auto array_index = 0; array_index < size; array_index++) {
            position_mid[array_index] = xgrid[array_index + 1];
            position_back[array_index] = xgrid[array_index];
            position_forward[array_index] = xgrid[array_index + 2];
            index_mid[array_index] = get_refractive_index(xgrid[array_index + 1], ygrid[idy], z + 0.5 * dz);
            index_back[array_index] = get_refractive_index(xgrid[array_index], ygrid[idy], z + 0.5 * dz);
            index_forward[array_index] = get_refractive_index(xgrid[array_index + 2], ygrid[idy], z + 0.5 * dz);

            rhs_slice[array_index] = rhs[array_index + 1][idy];
        }
        auto solution = OperatorSuite::solve_system(position_mid, position_back, position_forward, index_mid,
                                                    index_back, index_forward, reference_index, k0, preFactorLHS,
                                                    rhs_slice, pmlxPtr);


        //store into solution vector
        for (int idx = 1; idx < numx - 1; idx++) {
            half_step[idx][idy] = solution[idx - 1];
        }
    }


    // solve for next field value.
    multi_array<std::complex<double>, 2> full_step(extents[numx][numy]);
    std::fill_n(full_step.data(), full_step.num_elements(), 0.0);
    //invert (1+...Gy) operator, so only coupling between neighbours in y direction, this means we have numx independent problems
    for (int idx = 1; idx < numx - 1; idx++) {
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
            index_mid[array_index] = get_refractive_index(xgrid[idx], ygrid[array_index + 1], z + dz);
            index_back[array_index] = get_refractive_index(xgrid[idx], ygrid[array_index], z + dz);
            index_forward[array_index] = get_refractive_index(xgrid[idx], ygrid[array_index + 2], z + dz);
            //using the half_step now to fill the rhs vector.
            rhs_slice[array_index] = half_step[idx][array_index + 1];
        }

        auto solution = OperatorSuite::solve_system(position_mid, position_back, position_forward, index_mid,
                                                    index_back, index_forward, reference_index, k0, preFactorLHS,
                                                    rhs_slice, pmlyPtr);


        //store into solution vector
        for (int idy = 1; idy < numy - 1; idy++) {
            full_step[idx][idy] = solution[idy - 1];
        }
    }
    return full_step;
}


multi_array<double, 2> Solver::get_intensity(const multi_array<std::complex<double>, 2> &field) const {
    multi_array<double, 2> intensity(extents[numx][numy]);
    double max = 0.0;
    for (auto idx = 0; idx < field.shape()[0]; idx++) {
        for (auto idy = 0; idy < field.shape()[1]; idy++) {
            const double real = field[idx][idy].real();
            const double imag = field[idx][idy].imag();
            const double value = real * real + imag * imag;
            intensity[idx][idy] = value;
            max = std::max(max, value);
        }
    }
    for (auto idx = 0; idx < field.shape()[0]; idx++) {
        for (auto idy = 0; idy < field.shape()[1]; idy++) {
            intensity[idx][idy] /= max;
        }
    }
    return intensity;
}
