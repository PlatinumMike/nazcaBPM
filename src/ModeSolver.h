//
// Created by mike on 11/11/24.
//

#ifndef MODESOLVER_H
#define MODESOLVER_H

#include <filesystem>
#include "Port.h"
#include "Solver.h"


// specify how much output to write.
enum logging_level {
    ERROR, WARNING, INFO, DEBUG
};


class ModeSolver : public Solver {
public:
    ModeSolver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
               const Port &port,
               double scheme_parameter, double k0, double reference_index,
               const std::filesystem::path &absolute_path_output,
               logging_level level);

    /**
     * This will run the Solver to search for modes.
     * Currently only the fundamental mode is supported.
     * The iterations stop when a certain tolerance is reached, or when the maximum number of iterations is exceeded.
     * @param increment_x step size in x direction.
     * @param max_iterations maximum number of iterations allowed before the mode solver gives up.
     * This is typically the same as for the main BPM simulation, but you can specify a different step size.
     */
    void run(double increment_x, int max_iterations = 1000);

    /**
     * Get the mode overlap
     * @param bpm_field Field of the BPM simulation, already assumed to be interpolated on the grid of the ModeSolver.
     * @return coupling coefficient (mode overlap integral)
     */
    [[nodiscard]] double get_mode_overlap(const multi_array<std::complex<double>, 2> &bpm_field) const;

    /**
     * Get propagation constant of the mode that is found.
     * @return propagation constant (1/um).
     */
    [[nodiscard]] double get_beta() const;

    /**
      * Get effective index of the mode that is found.
      * @return effective index.
      */
    [[nodiscard]] double get_neff() const;

private:
    Port port;
    double beta;
    double neff;
    const double abs_tolerance = 1.0e-4;
    const int min_iterations = 10;
    const std::filesystem::path absolute_path_output;
    const logging_level level;

    /**
     * Rescales the field such that the integral of |u|^2 over the yz plane is 1.
     */
    void normalize_field(multi_array<std::complex<double>, 2> &field) const;

    double compute_beta(const multi_array<std::complex<double>, 2> &old_field,
                        const multi_array<std::complex<double>, 2> &new_field, double increment_x) const;

    double get_field_log(const multi_array<std::complex<double>, 2> &field) const;

    double get_norm(const multi_array<std::complex<double>, 2> &field) const;

    /**
        * Generates a 2D Gaussian profile
        * @param ygrid grid of y positions
        * @param zgrid grid of z positions
        * @param y0 center of Gaussian in y
        * @param z0 center of Gaussian in z
        * @param std_y standard deviation in y
        * @param std_z standard deviation in z
        * @return
        */
    multi_array<std::complex<double>, 2> get_initial_profile(const std::vector<double> &ygrid,
                                                             const std::vector<double> &zgrid, double y0,
                                                             double z0, double std_y,
                                                             double std_z) const;
};


#endif //MODESOLVER_H
