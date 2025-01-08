//
// Created by mike on 11/11/24.
//

#ifndef BPMSOLVER_H
#define BPMSOLVER_H

#include <filesystem>

#include "Solver.h"
#include "RectangularGrid3D.h"

class BpmSolver : public Solver {
public:
    BpmSolver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
              const RectangularGrid3D &grid,
              double scheme_parameter, double k0, double reference_index,
              const std::filesystem::path &absolute_path_output);

    /**
     * Run the BPM simulation.
     * It stops when the x position has reached the end, i.e. x=xmax.
     * @param initial_field Starting field at x=xmin.
     * @param field_slice_y y position to slice domain at
     * @param field_slice_z z position to slice domain at
     */
    void run(const multi_array<std::complex<double>, 2> &initial_field, const std::vector<double> &print_progress,
             double field_slice_y, double field_slice_z);

private:
    const RectangularGrid3D grid3d;
    const std::filesystem::path absolute_path_output;
};


#endif //BPMSOLVER_H
