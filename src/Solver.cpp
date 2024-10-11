//
// Created by mike on 10/11/24.
//

#include "Solver.h"
#include "AuxiliaryFunctions.h"
#include "IO/Readers.h"
#include "IO/hdf_writer.h"

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
               int numz, const double scheme_parameter): xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin),
                                                         zmax(zmax), numx(numx),
                                                         numy(numy), numz(numz), scheme_parameter(scheme_parameter) {
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
    // multi_array<double, 3> index_dataset(extents[numx][numy][numz]);
    // for (int i = 0; i < numx; i++) {
    //     for (int j = 0; j < numy; j++) {
    //         for (int k = 0; k < numz; k++) {
    //             index_dataset[i][j][k] = get_refractive_index(xgrid[i], ygrid[j], zgrid[k]);
    //         }
    //     }
    // }
    // auto xgrid_m = vector_to_multi_array(xgrid);
    // auto ygrid_m = vector_to_multi_array(ygrid);
    // auto zgrid_m = vector_to_multi_array(zgrid);
    //
    // // write to HDF5 format
    // H5::H5File file("index_data.h5", H5F_ACC_TRUNC);
    // write_hdf5(file, "refractive_index", index_dataset);
    // write_hdf5(file, "xgrid", xgrid_m);
    // write_hdf5(file, "ygrid", ygrid_m);
    // write_hdf5(file, "zgrid", zgrid_m);
    // file.close();

    multi_array<std::complex<double>, 2> field_slice_yz(extents[numy][numz]);
    record_slice(field, field_slice_yz, 0);


    const int index_1percent = numz / 100;
    const int index_10percent = numz / 10;
    const int index_50percent = numz / 2;

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int z_step = 1; z_step < numz; z_step++) {
        double current_z = zgrid[z_step - 1];
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


// multi_array<std::complex<double>, 2> Solver::get_derivative(const multi_array<std::complex<double>, 2> &field,
//                                                             const std::vector<double> &xgrid,
//                                                             const std::vector<double> &ygrid, double z) const {
//     multi_array<std::complex<double>, 2> derivative(extents[numx][numy]);
//     //update just the middle of the grid, because the boundary values will remain at zero. So set to zero once, then keep it there.
//     std::complex<double> prefactor = {0.0, -0.5 / beta_ref};
//     std::complex<double> qfactor_x_mid, qfactor_x_bottom, qfactor_x_top, qfactor_y_mid, qfactor_y_left, qfactor_y_right,
//             field_mid, field_left, field_right, field_top, field_bottom;
//     std::complex<double> value = {0.0, 0.0};
//     const double dx = xgrid[1] - xgrid[0];
//     const double dy = ygrid[1] - ygrid[0];
//     double index;
//     for (int idx = 0; idx < numx; idx++) {
//         for (int idy = 0; idy < numy; idy++) {
//             if (idx == 0 || idx == numx - 1 || idy == 0 || idy == numy - 1) {
//                 derivative[idx][idy] = std::complex<double>{0.0, 0.0};
//             } else {
//                 index = get_refractive_index(xgrid[idx], ygrid[idy], z);
//                 qfactor_x_mid = get_pml_factor1(xgrid[idx], ygrid[idy], z);
//                 qfactor_x_bottom = get_pml_factor1(xgrid[idx - 1], ygrid[idy], z);
//                 qfactor_x_top = get_pml_factor1(xgrid[idx + 1], ygrid[idy], z);
//                 qfactor_y_mid = get_pml_factor2(xgrid[idx], ygrid[idy], z);
//                 qfactor_y_left = get_pml_factor2(xgrid[idx], ygrid[idy - 1], z);
//                 qfactor_y_right = get_pml_factor2(xgrid[idx], ygrid[idy + 1], z);
//                 field_mid = field[idx][idy];
//                 field_bottom = field[idx - 1][idy];
//                 field_top = field[idx + 1][idy];
//                 field_left = field[idx][idy - 1];
//                 field_right = field[idx][idy + 1];
//
//                 //k0^2(n^2-n0^2)*u
//                 value = k0 * k0 * (index * index - reference_index * reference_index) * field_mid;
//                 //add Laplacian
//                 value += qfactor_x_mid / (2.0 * dx * dx) * (
//                     (qfactor_x_mid + qfactor_x_bottom) * field_bottom - (
//                         qfactor_x_top + 2.0 * qfactor_x_mid + qfactor_x_bottom) * field_mid + (
//                         qfactor_x_mid + qfactor_x_top) * field_top);
//                 value += qfactor_y_mid / (2.0 * dy * dy) * (
//                     (qfactor_y_mid + qfactor_y_left) * field_left - (
//                         qfactor_y_right + 2.0 * qfactor_y_mid + qfactor_y_left) * field_mid + (
//                         qfactor_y_mid + qfactor_y_right) * field_right);
//
//                 derivative[idx][idy] = prefactor * value;
//             }
//         }
//     }
//
//     return derivative;
// }

// multi_array<std::complex<double>, 2> Engine::do_step_rk4(const multi_array<std::complex<double>, 2> &field,
//                                                          const std::vector<double> &xgrid,
//                                                          const std::vector<double> &ygrid, double z,
//                                                          const double dz) const {
//     //The RK4 method can be rewritten to use less memory. Initialize the new_field as a copy of the old_field.
//     // Then compute k1, and add it to the new field. Also use it for k2, then deallocate k1. repeat for k2, k3, k4.
//     multi_array<std::complex<double>, 2> new_field(extents[numx][numy]);
//     multi_array<std::complex<double>, 2> temp_field(extents[numx][numy]);
//     const multi_array<std::complex<double>, 2> k1 = get_derivative(field, xgrid, ygrid, z);
//     for (int idx = 0; idx < numx; idx++) {
//         for (int idy = 0; idy < numy; idy++) {
//             temp_field[idx][idy] = field[idx][idy] + 0.5 * k1[idx][idy] * dz;
//         }
//     }
//     const multi_array<std::complex<double>, 2> k2 = get_derivative(temp_field, xgrid, ygrid, z + 0.5 * dz);
//     for (int idx = 0; idx < numx; idx++) {
//         for (int idy = 0; idy < numy; idy++) {
//             temp_field[idx][idy] = field[idx][idy] + 0.5 * k2[idx][idy] * dz;
//         }
//     }
//     const multi_array<std::complex<double>, 2> k3 = get_derivative(temp_field, xgrid, ygrid, z + 0.5 * dz);
//     for (int idx = 0; idx < numx; idx++) {
//         for (int idy = 0; idy < numy; idy++) {
//             temp_field[idx][idy] = field[idx][idy] + k3[idx][idy] * dz;
//         }
//     }
//     const multi_array<std::complex<double>, 2> k4 = get_derivative(temp_field, xgrid, ygrid, z + dz);
//
//
//     for (int idx = 0; idx < numx; idx++) {
//         for (int idy = 0; idy < numy; idy++) {
//             new_field[idx][idy] = field[idx][idy] + (
//                                       k1[idx][idy] + 2.0 * k2[idx][idy] + 2.0 * k3[idx][idy] + k4[idx][idy]) * dz / 6.0;
//         }
//     }
//
//
//     return new_field;
// }

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

multi_array<std::complex<double>, 2> Solver::do_step_cn(multi_array<std::complex<double>, 2> field,
    std::vector<double> xgrid, std::vector<double> ygrid, double z, double dz) {
    //todo
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
