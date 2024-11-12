//
// Created by mike on 11/11/24.
//

#include "ModeSolver.h"

#include <iostream>
#include <format>

ModeSolver::ModeSolver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source,
                       const RectangularGrid &grid, const double scheme_parameter, const double k0,
                       const double reference_index, const Port &port): Solver(
                                                                            geometry, pmlx, pmly, grid,
                                                                            scheme_parameter, k0, reference_index,
                                                                            true),
                                                                        port(port), sourcePtr(&source) {
    beta = 0.0;
    neff = 0.0;
}

void ModeSolver::run(const int max_iterations) {
    auto ygrid = gridPtr->get_ygrid();
    auto zgrid = gridPtr->get_zgrid();

    std::cout << std::format("Mode solving for port {} now...", port.get_name()) << std::endl;

    // slight offset so it excites also odd modes, this can be used to find the higher order modes
    constexpr double eps_y = 0.5;
    constexpr double eps_z = 0.3;
    //define field in current slice to be "field". Since it is a scalar BPM there is no polarization.
    internal_field = sourcePtr->get_initial_profile(
        ygrid, zgrid, port.get_y0() + eps_y, port.get_z0() + eps_z, 1.0, 1.0);

    double beta_old = beta;
    bool mode_found = false;
    double iterations_used = 0;

    for (int x_step = 0; x_step < max_iterations; x_step++) {
        auto new_field = do_step_cn(internal_field, port.get_x0(), gridPtr->get_dx());
        beta = compute_beta(internal_field, new_field);
        internal_field = new_field;
        if (x_step + 1 >= min_iterations && std::abs(beta - beta_old) < abs_tolerance) {
            //converged, so exiting the loop
            mode_found = true;
            iterations_used = x_step + 1;
        }
        beta_old = beta;

        // Normalize fields to have power=1. In principle we can just do this once at the end,
        // but in order to avoid numerical over/underflows of the Im-Dis method (exponential growth/decay), we do it every iteration.
        normalize_field(internal_field);

        if (mode_found) {
            break;
        }
    }
    neff = beta / k0;
    if (mode_found) {
        std::cout << std::format("Mode found after {} iterations, beta = {}, neff = {}", iterations_used, beta, neff) <<
                std::endl;
    } else {
        std::cout << std::format("Max iterations reached, mode not found, latest beta = {}, neff = {}", beta, neff) <<
                std::endl;
    }
}


double ModeSolver::get_field_log(const multi_array<std::complex<double>, 2> &field) const {
    std::complex<double> integral = {0, 0};
    for (auto i = 0; i < field.shape()[0]; i++) {
        for (auto j = 0; j < field.shape()[1]; j++) {
            integral += field[i][j];
        }
    }
    integral = std::log(integral);
    //the imaginary part should be negligible, so dropping it.
    return integral.real();
}

double ModeSolver::get_norm(const multi_array<std::complex<double>, 2> &field) const {
    // the norm is simply the integral of |u|^2 over y,z.
    // We approximate this as the sum of the grid values squared, times the surface are of a single cell.
    double norm = 0.0;
    for (auto i = 0; i < field.shape()[0]; i++) {
        for (auto j = 0; j < field.shape()[1]; j++) {
            norm += field[i][j].real() * field[i][j].real()
                    + field[i][j].imag() * field[i][j].imag();
        }
    }
    norm *= gridPtr->get_dy() * gridPtr->get_dz();
    return norm;
}

void ModeSolver::normalize_field(multi_array<std::complex<double>, 2> &field) const {
    const double norm = get_norm(field);
    for (auto i = 0; i < field.shape()[0]; i++) {
        for (auto j = 0; j < field.shape()[1]; j++) {
            field[i][j] /= norm;
        }
    }
}

double ModeSolver::get_mode_overlap(const multi_array<std::complex<double>, 2> &bpm_field) const {
    // Simplified formula for the case where we just have a single field component to work with.
    // this is justified for scalar and semi-vectorial BPM.
    const double norm1 = get_norm(internal_field);
    const double norm2 = get_norm(bpm_field);
    std::complex<double> cross_term = {0, 0};
    for (auto i = 0; i < internal_field.shape()[0]; i++) {
        for (auto j = 0; j < internal_field.shape()[1]; j++) {
            cross_term += std::conj(internal_field[i][j]) * bpm_field[i][j];
        }
    }
    cross_term *= gridPtr->get_dy() * gridPtr->get_dz();
    const double numerator = cross_term.real() * cross_term.real() + cross_term.imag() * cross_term.imag();
    return numerator / (norm1 * norm2);
}

double ModeSolver::compute_beta(const multi_array<std::complex<double>, 2> &old_field,
                                const multi_array<std::complex<double>, 2> &new_field) const {
    const double field_log_old = get_field_log(old_field);
    const double field_log_new = get_field_log(new_field);
    return reference_index * k0 + (field_log_new - field_log_old) / gridPtr->get_dx();
}

double ModeSolver::get_beta() const {
    return beta;
}

double ModeSolver::get_neff() const {
    return neff;
}
