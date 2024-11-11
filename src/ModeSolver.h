//
// Created by mike on 11/11/24.
//

#ifndef MODESOLVER_H
#define MODESOLVER_H
#include "Port.h"
#include "Solver.h"
#include "ModeHandler.h"


class ModeSolver : public Solver {
public:
    ModeSolver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source,
               const RectangularGrid &grid,
               double scheme_parameter, double k0, double reference_index, const Port &port);

    /**
     * This will run the Solver to search for modes.
     * Currently only the fundamental mode is supported.
     * The iterations stop when a certain tolerance is reached, or when the maximum number of iterations is exceeded.
     */
    void run(int max_iterations = 1000);

    //todo: this takes in a grid and field, and interpolates on the internal grid then does the mode overlap with the internal field.
    std::complex<double> get_mode_overlap(const multi_array<std::complex<double>, 2> &bpm_field) const;

    /**
     * Get propagation constant of the mode that is found.
     * @return propagation constant (1/um).
     */
    double get_beta() const;

    /**
     * Interpolate the internal field on a given new grid.
     * @param ygrid new y coordinates
     * @param zgrid new z coordinates
     * @return field on the new coordinates.
     */
    multi_array<std::complex<double>, 2> interpolate_field(const std::vector<double> &ygrid,
                                                           const std::vector<double> &zgrid) const;

private:
    Port port;
    const ModeHandler *sourcePtr;
    double beta;

    multi_array<std::complex<double>, 2> internal_field;
};


#endif //MODESOLVER_H
