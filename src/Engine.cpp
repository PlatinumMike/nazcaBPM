//
// Created by mike on 9/29/24.
//

#include "Engine.h"
#include "AuxiliaryFunctions.h"
#include "Writers.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <chrono>


Engine::Engine() {
    std::cout << "Launching new engine now" << std::endl;
}

void Engine::run() {
    numx = static_cast<int>(domain_len_x * resolution_x);
    numy = static_cast<int>(domain_len_y * resolution_y);
    // int num_slice = numx * numy;
    double xgrid[numx];
    double ygrid[numy];
    //TODO: replace this with a function that returns a vector, that is a bit safer than a raw array.
    AuxiliaryFunctions::linspace(xgrid, -0.5 * domain_len_x, 0.5 * domain_len_x, numx);
    AuxiliaryFunctions::linspace(ygrid, -0.5 * domain_len_y, 0.5 * domain_len_y, numy);
    const double dx = xgrid[1] - xgrid[0];
    const double dy = ygrid[1] - ygrid[0];
    double constexpr safety_factor = 0.5; // best to use a slightly smaller dz step to be sure it is stable.
    numz = static_cast<int>(domain_len_z / (safety_factor * get_min_dz(dx, dy)));
    double zgrid[numz];
    AuxiliaryFunctions::linspace(zgrid, 0.0, domain_len_y, numz);
    const double dz = zgrid[1] - zgrid[0];

    std::cout << "Lx = " << domain_len_x << ", Ly = " << domain_len_y << ", Lz = " << domain_len_z << std::endl;
    std::cout << "Nx = " << numx << ", Ny = " << numy << ", Nz = " << numz << std::endl;
    std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;

    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    std::vector<std::complex<double> > field = get_initial_profile(xgrid, ygrid);

    Writers::write_to_csv(field, "initial_field.csv", numx, numy);

    std::vector<std::complex<double> > field_slice_yz;
    record_slice(field, field_slice_yz);


    int index_1percent = numz / 100;
    int index_10percent = numz / 10;
    int index_50percent = numz / 2;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int z_step = 1; z_step < numz; z_step++) {
        double current_z = zgrid[z_step - 1];
        field = do_step_rk4(field, xgrid, ygrid, current_z, dz);
        record_slice(field, field_slice_yz);
        //printing rough indication of simulation progress
        if (z_step == index_1percent) {
            std::cout << "1% reached" << std::endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<
                    " (s)" <<
                    std::endl;
        }
        if (z_step == index_10percent) {
            std::cout << "10% reached" << std::endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<
                    " (s)" <<
                    std::endl;
        }
        if (z_step == index_50percent) {
            std::cout << "50% reached" << std::endl;
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<
                    " (s)" <<
                    std::endl;
        }
    }
    std::cout << "Engine run completed" << std::endl;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " (s)" <<
            std::endl;

    Writers::write_to_csv(field, "final_field.csv", numx, numy);
    Writers::write_to_csv(field_slice_yz, "field_slice.csv", numz, numy);
}


double Engine::get_refractive_index(const double x, const double y, double z) const {
    //TODO: replace this with proper list of polygons to search in. Nazca export.
    double wg_width = 1.000;
    double wg_height = 0.350;
    double x0 = 0.0;
    double y0 = 0.0;
    if (x < x0 + 0.5 * wg_height && x > x0 - 0.5 * wg_height && y < y0 + 0.5 * wg_width && y > y0 - 0.5 * wg_width) {
        return core_index;
    } else {
        return background_index;
    }
}

double Engine::get_min_dz(double dx, double dy) const {
    // from this paper "Analysis of 2-Invariant and 2-Variant Semiconductor Rib Waveguides by Explicit Finite Difference BeamPropagation Method with Nonuniform MeshConfiguration", 1991.
    double max_index = core_index; // max index in the whole domain.
    double denom = 4.0 / (dx * dx) + 4.0 / (dy * dy) + k0 * k0 * (
                       max_index * max_index - reference_index * reference_index);
    double dz_min = 2.0 * beta_ref / denom;
    return dz_min;
}

// just some mock-up profile. I should use the actual fundamental mode.
std::vector<std::complex<double> > Engine::get_initial_profile(const double *xgrid, const double *ygrid) const {
    //todo: in the future this would be read from a file, which is generated by a mode solver, e.g. WGMS3D.
    std::vector<std::complex<double> > initial_profile;
    std::complex<double> value = {0.0, 0.0};
    double x0 = 0.0;
    double y0 = 0.0;
    double x = 0.0;
    double y = 0.0;
    double std_x = 1.0; //standard deviation
    double std_y = 1.0;
    double deltax, deltay;
    for (int idx = 0; idx < numx; idx++) {
        for (int idy = 0; idy < numy; idy++) {
            if (idx == 0 || idx == numx - 1 || idy == 0 || idy == numy - 1) {
                //set boundary values to zero. Simple Dirichlet BCs. Reflections are negligible since a PML will be placed in front of the metal wall anyway.
                value = {0.0, 0.0};
            } else {
                x = xgrid[idx];
                y = ygrid[idy];
                deltax = (x - x0) / std_x;
                deltay = (y - y0) / std_y;
                value = {std::exp(-0.5 * deltax * deltax - 0.5 * deltay * deltay), 0.0};
            }
            initial_profile.push_back(value);
        }
    }
    return initial_profile;
}


std::vector<std::complex<double> > Engine::get_derivative(const std::vector<std::complex<double> > &field,
                                                          double *xgrid,
                                                          double *ygrid, double z) const {
    std::vector<std::complex<double> > derivative;
    //update just the middle of the grid, because the boundary values will remain at zero. So set to zero once, then keep it there.
    std::complex<double> prefactor = {0.0, -0.5 / beta_ref};
    std::complex<double> qfactor_x_mid, qfactor_x_bottom, qfactor_x_top, qfactor_y_mid, qfactor_y_left, qfactor_y_right,
            field_mid, field_left, field_right, field_top, field_bottom;
    std::complex<double> value = {0.0, 0.0};
    double dx = xgrid[1] - xgrid[0];
    double dy = ygrid[1] - ygrid[0];
    double index;
    for (int idx = 0; idx < numx; idx++) {
        for (int idy = 0; idy < numy; idy++) {
            if (idx == 0 || idx == numx - 1 || idy == 0 || idy == numy - 1) {
                derivative.push_back(std::complex<double>{0.0, 0.0});
            } else {
                //TODO: clean up this part, the selecting of indices is very error-prone.
                int i_mid = idx * numy + idy;
                int i_right = get_neighbour_index(idx, idy, numy, 0);
                int i_top = get_neighbour_index(idx, idy, numy, 1);
                int i_left = get_neighbour_index(idx, idy, numy, 2);
                int i_bottom = get_neighbour_index(idx, idy, numy, 3);
                index = get_refractive_index(xgrid[idx], ygrid[idy], z);
                qfactor_x_mid = get_qfactorx(xgrid[idx], ygrid[idy], z);
                qfactor_x_bottom = get_qfactorx(xgrid[idx - 1], ygrid[idy], z);
                qfactor_x_top = get_qfactorx(xgrid[idx + 1], ygrid[idy], z);
                qfactor_y_mid = get_qfactory(xgrid[idx], ygrid[idy], z);
                qfactor_y_left = get_qfactory(xgrid[idx], ygrid[idy - 1], z);
                qfactor_y_right = get_qfactory(xgrid[idx], ygrid[idy + 1], z);
                field_mid = field[i_mid];
                field_left = field[i_left];
                field_right = field[i_right];
                field_top = field[i_top];
                field_bottom = field[i_bottom];
                //k0^2(n^2-n0^2)*u
                value = k0 * k0 * (index * index - reference_index * reference_index) * field_mid;
                //add Laplacian
                value += qfactor_x_mid / (2.0 * dx * dx) * (
                    (qfactor_x_mid + qfactor_x_bottom) * field_bottom - (
                        qfactor_x_top + 2.0 * qfactor_x_mid + qfactor_x_bottom) * field_mid + (
                        qfactor_x_mid + qfactor_x_top) * field_top);
                value += qfactor_y_mid / (2.0 * dy * dy) * (
                    (qfactor_y_mid + qfactor_y_left) * field_left - (
                        qfactor_y_right + 2.0 * qfactor_y_mid + qfactor_y_left) * field_mid + (
                        qfactor_y_mid + qfactor_y_right) * field_right);

                derivative.push_back(prefactor * value);
            }
        }
    }

    return derivative;
}

std::vector<std::complex<double> > Engine::do_step_euler(const std::vector<std::complex<double> > &field, double *xgrid,
                                                         double *ygrid, double z, double dz) const {
    std::vector<std::complex<double> > new_field;
    std::vector<std::complex<double> > k1 = get_derivative(field, xgrid, ygrid, z);
    for (size_t index = 0; index < field.size(); index++) {
        new_field.push_back(field[index] + k1[index] * dz);
    }

    return new_field;
}

std::vector<std::complex<double> > Engine::do_step_rk4(const std::vector<std::complex<double> > &field, double *xgrid,
                                                       double *ygrid, double z, double dz) const {
    //The RK4 method can be rewritten to use less memory. Initialize the new_field as a copy of the old_field.
    // Then compute k1, and add it to the new field. Also use it for k2, then deallocate k1. repeat for k2, k3, k4.
    std::vector<std::complex<double> > new_field;
    std::vector<std::complex<double> > temp_field(field.size());
    std::vector<std::complex<double> > k1 = get_derivative(field, xgrid, ygrid, z);
    for (size_t index = 0; index < field.size(); index++) {
        temp_field[index] = field[index] + 0.5 * k1[index] * dz;
    }
    std::vector<std::complex<double> > k2 = get_derivative(temp_field, xgrid, ygrid, z + 0.5 * dz);
    for (size_t index = 0; index < field.size(); index++) {
        temp_field[index] = field[index] + 0.5 * k2[index] * dz;
    }
    std::vector<std::complex<double> > k3 = get_derivative(temp_field, xgrid, ygrid, z + 0.5 * dz);
    for (size_t index = 0; index < field.size(); index++) {
        temp_field[index] = field[index] + k3[index] * dz;
    }
    std::vector<std::complex<double> > k4 = get_derivative(temp_field, xgrid, ygrid, z + dz);


    for (size_t index = 0; index < field.size(); index++) {
        new_field.push_back(field[index] + (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index]) * dz / 6.0);
    }

    return new_field;
}


double Engine::get_conductivity_base(double x, double xmin, double xmax) const {
    //second order polynomial PML strength
    if (x < xmin + pml_thickness) {
        double base_value = ((xmin + pml_thickness) - x) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else if (x > xmax - pml_thickness) {
        double base_value = (x - (xmax - pml_thickness)) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else {
        return 0.0;
    }
}

int Engine::get_neighbour_index(int i, int j, int ny, int direction) const {
    // direction = 0 means to the right, so increment the y coordinate.
    // direction = 1 means to the top, so increment the x coordinate.
    // direction = 2 means to the left, so decrement the y coordinate.
    // direction = 3 means to the bottom, so decrement the x coordinate.
    int base_index = i * ny + j;
    if (direction == 0) {
        return base_index + 1;
    } else if (direction == 1) {
        return base_index + ny;
    } else if (direction == 2) {
        return base_index - 1;
    } else {
        return base_index - ny;
    }
}

void Engine::record_slice(const std::vector<std::complex<double> > &buffer,
                          std::vector<std::complex<double> > &storage) {
    int index_xmid = numx / 2; //yz slice at x=0 approximately
    for (int idy = 0; idy < numy; idy++) {
        storage.push_back(buffer[index_xmid * numy + idy]);
    }

    // for (auto value : buffer) {
    //     storage.push_back(value);
    // }
}

double Engine::get_conductivityx(double x) const {
    return get_conductivity_base(x, -0.5 * domain_len_x, 0.5 * domain_len_x);
}

double Engine::get_conductivityy(double y) const {
    return get_conductivity_base(y, -0.5 * domain_len_y, 0.5 * domain_len_y);
}


std::complex<double> Engine::get_qfactorx(double x, double y, double z) const {
    double pml_index = get_refractive_index(x, y, z); // equal to bordering index value.
    return 1.0 / std::complex{1.0, -get_conductivityx(x) / (pml_index * pml_index)};
}

std::complex<double> Engine::get_qfactory(double x, double y, double z) const {
    double pml_index = get_refractive_index(x, y, z); // equal to bordering index value.
    return 1.0 / std::complex{1.0, -get_conductivityy(x) / (pml_index * pml_index)};
}

std::vector<double> Engine::get_intensity(const std::vector<std::complex<double> > &field) {
    std::vector<double> intensity;
    double real, imag;
    for (auto value: field) {
        real = value.real();
        imag = value.imag();
        intensity.push_back(real * real + imag * imag);
    }
    double max = *std::max_element(intensity.begin(), intensity.end());
    for (double &intens: intensity) {
        intens /= max;
    }
    return intensity;
}
