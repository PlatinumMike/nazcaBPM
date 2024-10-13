//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "AuxiliaryFunctions.h"
#include "IO/Readers.h"
#include "IO/hdf_writer.h"
#include "TriDiag.h"

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

    dump_index_slice("index_start.h5", 'z', 0.0, xgrid, ygrid);
    dump_index_slice("index_end.h5", 'z', zgrid[numz - 1], xgrid, ygrid);
    dump_index_slice("index_cross.h5", 'x', 0.0, ygrid, zgrid);

    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    multi_array<std::complex<double>, 2> field = sourcePtr->get_initial_profile(xgrid, ygrid, 0.0, 0.0, 1.0, 1.0);

    write_cmplx_hdf5("initial_field.h5", field);

    //duming index on entire 3D grid into hdf5 using BOOST
    multi_array<double, 3> index_dataset(extents[numx][numy][numz]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            for (int k = 0; k < numz; k++) {
                index_dataset[i][j][k] = get_refractive_index(xgrid[i], ygrid[j], zgrid[k]);
            }
        }
    }
    auto xgrid_m = vector_to_multi_array(xgrid);
    auto ygrid_m = vector_to_multi_array(ygrid);
    auto zgrid_m = vector_to_multi_array(zgrid);

    // write to HDF5 format
    H5::H5File file("index_data.h5", H5F_ACC_TRUNC);
    write_hdf5(file, "refractive_index", index_dataset);
    write_hdf5(file, "xgrid", xgrid_m);
    write_hdf5(file, "ygrid", ygrid_m);
    write_hdf5(file, "zgrid", zgrid_m);
    file.close();

    multi_array<std::complex<double>, 2> field_slice_yz(extents[numy][numz]);
    record_slice(field, field_slice_yz, 0);


    const int index_1percent = numz / 100;
    const int index_10percent = numz / 10;
    const int index_50percent = numz / 2;

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int z_step = 1; z_step < numz; z_step++) {
        const double current_z = zgrid[z_step - 1];
        field = do_step_cn(field, xgrid, ygrid, current_z, dz);
        record_slice(field, field_slice_yz, z_step);
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

    write_cmplx_hdf5("final_field.h5", field);
    write_cmplx_hdf5("field_slice.h5", field_slice_yz);
}

double Solver::get_refractive_index(const double x_bpm, const double y_bpm, double z_bpm) const {
    //convert to nazca coordinates
    return geometryPtr->get_index(z_bpm, -y_bpm, x_bpm);
}


void Solver::record_slice(const multi_array<std::complex<double>, 2> &buffer,
                          multi_array<std::complex<double>, 2> &storage, int idz) const {
    const int index_xmid = numx / 2; //yz slice at x=0 approximately
    for (int idy = 0; idy < numy; idy++) {
        storage[idy][idz] = buffer[index_xmid][idy];
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

multi_array<double, 1> Solver::vector_to_multi_array(const std::vector<double> &vec) const {
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
                                                        const double dz) {
    const auto pre_factor = scheme_parameter * dz / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //get RHS vector, store as multi-array.
    auto rhs = get_rhs(field, xgrid, ygrid, z, dz);

    //solver for half step field
    multi_array<std::complex<double>, 2> half_step(extents[numx][numy]);
    std::fill_n(half_step.data(), half_step.num_elements(), 0.0);
    //invert (1+...Gx) operator, so only coupling between neighbours in x direction, this means we have numy independent problems
    for (int idy = 1; idy < numy - 1; idy++) {
        auto size = numx - 2;
        std::vector<double> xmid(size);
        std::vector<double> xneighbor_previous(size);
        std::vector<double> xneighbor_next(size);
        std::vector<double> ymid(size);
        std::vector<double> yneighbor_previous(size);
        std::vector<double> yneighbor_next(size);
        std::vector<std::complex<double> > rhs_slice(size);
        //fill vectors
        for (auto array_index = 0; array_index < size; array_index++) {
            xmid[array_index] = xgrid[array_index + 1];
            xneighbor_previous[array_index] = xgrid[array_index];
            xneighbor_next[array_index] = xgrid[array_index + 2];
            //derivative in x direction, so neighbors have the same y value
            ymid[array_index] = ygrid[idy];
            yneighbor_previous[array_index] = ygrid[idy];
            yneighbor_next[array_index] = ygrid[idy];
            rhs_slice[array_index] = rhs[array_index + 1][idy];
        }
        auto solution = solve_system(xmid, xneighbor_previous, xneighbor_next, ymid, yneighbor_previous, yneighbor_next,
                                     dx, z + 0.5 * dz, pre_factor, rhs_slice, pmlxPtr);


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
        std::vector<double> xmid(size);
        std::vector<double> xneighbor_previous(size);
        std::vector<double> xneighbor_next(size);
        std::vector<double> ymid(size);
        std::vector<double> yneighbor_previous(size);
        std::vector<double> yneighbor_next(size);
        std::vector<std::complex<double> > rhs_slice(size);
        //fill vectors
        for (auto array_index = 0; array_index < size; array_index++) {
            xmid[array_index] = xgrid[idx];
            xneighbor_previous[array_index] = xgrid[idx];
            xneighbor_next[array_index] = xgrid[idx];
            //derivative in x direction, so neighbors have the same y value
            ymid[array_index] = ygrid[array_index + 1];
            yneighbor_previous[array_index] = ygrid[array_index];
            yneighbor_next[array_index] = ygrid[array_index + 2];
            //using the half_step now to fill the rhs vector.
            rhs_slice[array_index] = half_step[idx][array_index + 1];
        }
        auto solution = solve_system(xmid, xneighbor_previous, xneighbor_next, ymid, yneighbor_previous, yneighbor_next,
                                     dy, z + dz, pre_factor, rhs_slice, pmlyPtr);


        //store into solution vector
        for (int idy = 1; idy < numy - 1; idy++) {
            full_step[idx][idy] = solution[idy - 1];
        }
    }
    return full_step;
}

multi_array<std::complex<double>, 2> Solver::get_rhs(const multi_array<std::complex<double>, 2> &field,
                                                     const std::vector<double> &xgrid, const std::vector<double> &ygrid,
                                                     const double z, const double dz) const {
    // prefactor is (alpha-1)*dz*i/(2k0*n0).
    const std::complex<double> pre_factor =
            (scheme_parameter - 1.0) * dz / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};
    //right hand side vector
    multi_array<std::complex<double>, 2> rhs1(extents[numx][numy]);
    // apply the operator (1-...Gy)
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                rhs1[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // x value on the mid-point and neighbors is the same, as we only apply the derivative in y direction here.
                rhs1[i][j] = apply_right_hand_operator(xgrid[i], xgrid[i], xgrid[i], ygrid[j], ygrid[j - 1],
                                                       ygrid[j + 1],
                                                       dy, z,
                                                       pre_factor, field[i][j], field[i][j - 1], field[i][j + 1],
                                                       pmlyPtr);
            }
        }
    }
    //now apply the operator (1-...Gx)
    multi_array<std::complex<double>, 2> rhs2(extents[numx][numy]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                rhs2[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // y value on the mid-point and neighbors is the same, as we only apply the derivative in x direction here.
                rhs2[i][j] = apply_right_hand_operator(xgrid[i], xgrid[i - 1], xgrid[i + 1], ygrid[j], ygrid[j],
                                                       ygrid[j],
                                                       dx, z,
                                                       pre_factor, rhs1[i][j], rhs1[i - 1][j], rhs1[i][j], pmlxPtr);
            }
        }
    }
    return rhs2;
}

std::complex<double> Solver::apply_right_hand_operator(const double xmid, const double xneighbor_previous,
                                                       const double xneighbor_next,
                                                       const double ymid,
                                                       const double yneighbor_previous, const double yneighbor_next,
                                                       const double neighbor_distance, const double z,
                                                       const std::complex<double> preFactor,
                                                       const std::complex<double> fieldValueMid,
                                                       const std::complex<double> fieldValue_previous,
                                                       const std::complex<double> fieldValue_next,
                                                       const PML *pmlPtr) const {
    const double index_mid = get_refractive_index(xmid, ymid, z);
    const double index_previous = get_refractive_index(xneighbor_previous, yneighbor_previous, z);
    const double index_next = get_refractive_index(xneighbor_next, yneighbor_next, z);
    auto pmlfactor_mid = pmlPtr->get_pml_factor(xmid, ymid, z, index_mid);
    auto pmlfactor_previous = pmlPtr->get_pml_factor(xneighbor_previous, yneighbor_previous, z, index_previous);
    auto pmlfactor_next = pmlPtr->get_pml_factor(xneighbor_next, yneighbor_next, z, index_next);

    std::complex<double> value = fieldValueMid;
    const std::complex<double> coefficient = preFactor * pmlfactor_mid / (2.0 * neighbor_distance * neighbor_distance);

    // add index contribution
    value += preFactor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index) * fieldValueMid;
    // add double derivative contribution
    value += coefficient * ((pmlfactor_mid + pmlfactor_next) * fieldValue_next - (
                                pmlfactor_next + 2.0 * pmlfactor_mid + pmlfactor_previous) * fieldValueMid + (
                                pmlfactor_mid + pmlfactor_previous) * fieldValue_previous);
    return value;
}

std::vector<std::complex<double> > Solver::solve_system(const std::vector<double> &xmid,
                                                        const std::vector<double> &xneighbor_previous,
                                                        const std::vector<double> &xneighbor_next,
                                                        const std::vector<double> &ymid,
                                                        const std::vector<double> &yneighbor_previous,
                                                        const std::vector<double> &yneighbor_next,
                                                        const double neighbor_distance, const double z,
                                                        const std::complex<double> preFactor,
                                                        const std::vector<std::complex<double> > &rhs_slice,
                                                        const PML *pmlPtr) const {
    const auto size = rhs_slice.size();
    // compute matrix entries
    std::vector<std::complex<double> > lower_diag(size - 3);
    std::vector<std::complex<double> > diag(size - 2);
    std::vector<std::complex<double> > upper_diag(size - 3);
    std::vector<std::complex<double> > temp_rhs(size - 2);
    for (auto idx = 1; idx < size - 1; idx++) {
        const double index_mid = get_refractive_index(xmid[idx], ymid[idx], z);
        const double index_previous = get_refractive_index(xneighbor_previous[idx], yneighbor_previous[idx], z);
        const double index_next = get_refractive_index(xneighbor_next[idx], yneighbor_next[idx], z);
        auto pmlfactor_mid = pmlPtr->get_pml_factor(xmid[idx], ymid[idx], z, index_mid);
        auto pmlfactor_previous = pmlPtr->get_pml_factor(xneighbor_previous[idx], yneighbor_previous[idx], z,
                                                         index_previous);
        auto pmlfactor_next = pmlPtr->get_pml_factor(xneighbor_next[idx], yneighbor_next[idx], z, index_next);
        const std::complex<double> coefficient =
                preFactor * pmlfactor_mid / (2.0 * neighbor_distance * neighbor_distance);


        diag[idx - 1] = 1.0 + preFactor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index)
                        - coefficient * (pmlfactor_next + 2.0 * pmlfactor_mid + pmlfactor_previous);
        temp_rhs[idx - 1] = rhs_slice[idx];
        if (idx > 1) {
            lower_diag[idx - 2] = coefficient * (pmlfactor_mid + pmlfactor_previous);
        }
        if (idx < size - 2) {
            upper_diag[idx - 1] = coefficient * (pmlfactor_mid + pmlfactor_next);
        }
    }
    //solve
    auto solution = TriDiag<std::complex<double> >::solve_thomas(lower_diag, diag, upper_diag, temp_rhs);

    return solution;
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
