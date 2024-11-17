//
// Created by mike on 11/11/24.
//

#ifndef BPMSOLVER_H
#define BPMSOLVER_H

#include "Solver.h"

class BpmSolver : public Solver {
public:
    BpmSolver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
              const RectangularGrid &grid,
              double scheme_parameter, double k0, double reference_index);

    /**
     * Run the BPM simulation.
     * It stops when the x position has reached the end, i.e. x=xmax.
     * @param initial_field Starting field at x=xmin.
     */
    void run(const multi_array<std::complex<double>, 2> &initial_field);
};


#endif //BPMSOLVER_H
