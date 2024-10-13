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
// todo: also, a lot of code duplicatation here in do_step_cn, and get_rhs. Generalize, then clean up.
multi_array<std::complex<double>, 2> Solver::do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                        const std::vector<double> &xgrid,
                                                        const std::vector<double> &ygrid, const double z,
                                                        const double dz) {
    const auto pre_factor = scheme_parameter * dz / (2.0 * k0 * reference_index) * std::complex<double>{0.0, 1.0};

    //get RHS vector, store as multi-array.
    auto rhs = get_rhs(field, xgrid, ygrid, z, dz);

    //solver for half step field
    multi_array<std::complex<double>, 2> half_step(extents[numx][numy]);
    std::fill(half_step.data(), half_step.data() + half_step.num_elements(), 0.0);
    //invert (1+...Gx) operator, so only coupling between neighbours in x direction, this means we have numy independent problems
    for (int idy = 1; idy < numy - 1; idy++) {
        // compute matrix entries
        std::vector<std::complex<double> > lower_diag;
        std::vector<std::complex<double> > diag;
        std::vector<std::complex<double> > upper_diag;
        std::vector<std::complex<double> > temp_rhs;
        for (int idx = 1; idx < numx - 1; idx++) {
            std::complex<double> entry_d = {1.0, 0.0};
            std::complex<double> entry_ld = {0.0, 0.0};
            std::complex<double> entry_ud = {0.0, 0.0};
            const double y = ygrid[idy];
            const double index_mid = get_refractive_index(xgrid[idx], y, z);
            const double index_previous = get_refractive_index(xgrid[idx - 1], y, z);
            const double index_next = get_refractive_index(xgrid[idx + 1], y, z);
            auto pmlfactor_x_mid = pmlxPtr->get_pml_factor(xgrid[idx], y, z, index_mid);
            auto pmlfactor_x_previous = pmlxPtr->get_pml_factor(xgrid[idx - 1], y, z, index_previous);
            auto pmlfactor_x_next = pmlxPtr->get_pml_factor(xgrid[idx + 1], y, z, index_next);

            entry_d += pre_factor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index);
            entry_d += -pre_factor * pmlfactor_x_mid / (2.0 * dx * dx) * (
                pmlfactor_x_next + 2.0 * pmlfactor_x_mid + pmlfactor_x_previous);
            entry_ld += pre_factor * pmlfactor_x_mid / (2.0 * dx * dx) * (pmlfactor_x_mid + pmlfactor_x_previous);
            entry_ud += pre_factor * pmlfactor_x_mid / (2.0 * dx * dx) * (pmlfactor_x_mid + pmlfactor_x_next);

            diag.push_back(entry_d);
            temp_rhs.push_back(rhs[idx][idy]);
            if (idx > 1) {
                lower_diag.push_back(entry_ld);
            }
            if (idx < numx - 2) {
                upper_diag.push_back(entry_ud);
            }
        }
        //solve
        auto solution = TriDiag<std::complex<double> >::solve_thomas(lower_diag, diag, upper_diag, temp_rhs);
        //store into solution vector
        for (int idx = 1; idx < numx - 1; idx++) {
            half_step[idx][idy] = solution[idx];
        }
    }


    // solve for next field value.
    multi_array<std::complex<double>, 2> full_step(extents[numx][numy]);
    std::fill(full_step.data(), full_step.data() + full_step.num_elements(), 0.0);
    //invert (1+...Gy) operator, so only coupling between neighbours in y direction, this means we have numx independent problems
    for (int idx = 1; idx < numx - 1; idx++) {
        // compute matrix entries
        std::vector<std::complex<double> > lower_diag;
        std::vector<std::complex<double> > diag;
        std::vector<std::complex<double> > upper_diag;
        std::vector<std::complex<double> > temp_rhs;
        for (int idy = 1; idy < numy - 1; idy++) {
            std::complex<double> entry_d = {1.0, 0.0};
            std::complex<double> entry_ld = {0.0, 0.0};
            std::complex<double> entry_ud = {0.0, 0.0};
            const double x = xgrid[idx];
            const double index_mid = get_refractive_index(x, ygrid[idy], z);
            const double index_previous = get_refractive_index(x, ygrid[idy - 1], z);
            const double index_next = get_refractive_index(x, ygrid[idy + 1], z);
            auto pmlfactor_y_mid = pmlxPtr->get_pml_factor(x, ygrid[idy], z, index_mid);
            auto pmlfactor_y_previous = pmlxPtr->get_pml_factor(x, ygrid[idy - 1], z, index_previous);
            auto pmlfactor_y_next = pmlxPtr->get_pml_factor(x, ygrid[idy + 1], z, index_next);

            entry_d += pre_factor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index);
            entry_d += -pre_factor * pmlfactor_y_mid / (2.0 * dy * dy) * (
                pmlfactor_y_next + 2.0 * pmlfactor_y_next + pmlfactor_y_previous);
            entry_ld += pre_factor * pmlfactor_y_mid / (2.0 * dy * dy) * (pmlfactor_y_mid + pmlfactor_y_previous);
            entry_ud += pre_factor * pmlfactor_y_mid / (2.0 * dy * dy) * (pmlfactor_y_mid + pmlfactor_y_next);

            diag.push_back(entry_d);
            temp_rhs.push_back(half_step[idx][idy]);
            if (idy > 1) {
                lower_diag.push_back(entry_ld);
            }
            if (idy < numy - 2) {
                upper_diag.push_back(entry_ud);
            }
        }
        //solve
        auto solution = TriDiag<std::complex<double> >::solve_thomas(lower_diag, diag, upper_diag, temp_rhs);
        //store into solution vector
        for (int idy = 1; idy < numy - 1; idy++) {
            full_step[idx][idy] = solution[idx];
        }
    }
    return full_step;
}

multi_array<std::complex<double>, 2> Solver::get_rhs(const multi_array<std::complex<double>, 2> &field,
                                                     const std::vector<double> &xgrid, const std::vector<double> &ygrid,
                                                     const double z, const double dz) const {
    const std::complex<double> pre_factor = (scheme_parameter - 1.0) * dz / (2.0 * k0 * reference_index) * std::complex<
                                                double>
                                            {0.0, 1.0};
    //right hand side vector, make deep copy
    multi_array<std::complex<double>, 2> rhs1 = field;
    // apply the operator (1-...Gy)
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                rhs1[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // add Gy operator.
                const double x = xgrid[i];
                const double index_mid = get_refractive_index(x, ygrid[j], z);
                const double index_previous = get_refractive_index(x, ygrid[j - 1], z);
                const double index_next = get_refractive_index(x, ygrid[j + 1], z);
                auto pmlfactor_y_mid = pmlyPtr->get_pml_factor(x, ygrid[j], z, index_mid);
                auto pmlfactor_y_previous = pmlyPtr->get_pml_factor(x, ygrid[j - 1], z, index_previous);
                auto pmlfactor_y_next = pmlyPtr->get_pml_factor(x, ygrid[j + 1], z, index_next);
                // add index contribution
                rhs1[i][j] += pre_factor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index) *
                        field[i][j];
                // add double derivative contribution
                rhs1[i][j] += pre_factor * pmlfactor_y_mid / (2.0 * dy * dy) * (
                    (pmlfactor_y_mid + pmlfactor_y_next) * field[i][j + 1] - (
                        pmlfactor_y_next + 2.0 * pmlfactor_y_mid + pmlfactor_y_previous) * field[i][j] + (
                        pmlfactor_y_mid + pmlfactor_y_previous) * field[i][j - 1]);
            }
        }
    }
    //now apply the operator (1-...Gx)
    multi_array<std::complex<double>, 2> rhs2 = rhs1;
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                rhs2[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // add Gx operator.
                const double y = ygrid[j];
                const double index_mid = get_refractive_index(xgrid[i], y, z);
                const double index_previous = get_refractive_index(xgrid[i - 1], y, z);
                const double index_next = get_refractive_index(xgrid[i + 1], y, z);
                auto pmlfactor_x_mid = pmlxPtr->get_pml_factor(xgrid[i], y, z, index_mid);
                auto pmlfactor_x_previous = pmlxPtr->get_pml_factor(xgrid[i - 1], y, z, index_previous);
                auto pmlfactor_x_next = pmlxPtr->get_pml_factor(xgrid[i + 1], y, z, index_next);
                // add index contribution
                rhs2[i][j] += pre_factor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index) *
                        rhs1[
                            i][j];
                // add double derivative contribution
                rhs2[i][j] += pre_factor * pmlfactor_x_mid / (2.0 * dx * dx) * (
                    (pmlfactor_x_mid + pmlfactor_x_next) * rhs1[i + 1][j] - (
                        pmlfactor_x_next + 2.0 * pmlfactor_x_mid + pmlfactor_x_previous) * rhs1[i][j] + (
                        pmlfactor_x_mid + pmlfactor_x_previous) * rhs1[i - 1][j]);
            }
        }
    }
    return rhs2;
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
