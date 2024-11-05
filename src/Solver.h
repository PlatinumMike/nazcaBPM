//
// Created by mike on 10/11/24.
//

#ifndef SOLVER_H
#define SOLVER_H

#include <complex>
#include "Geometry.h"
#include <boost/multi_array.hpp>

#include "ModeHandler.h"
#include "PML.h"
#include "RectangularGrid.h"
using boost::multi_array;

// todo: for the current implementation many arrays are copied at every iteration.
// Performance can be improved by overwriting them in place, at least where that is safe.


class Solver {
public:
    Solver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source,
           const RectangularGrid &grid,
           double scheme_parameter, double k0, double reference_index);

    void run();

    multi_array<double, 2> get_intensity(const multi_array<std::complex<double>, 2> &field) const;

private:
    double scheme_parameter = 0.5;
    double k0;
    double reference_index;
    // pointers to simulation objects
    const Geometry *geometryPtr;
    const PML *pmlyPtr;
    const PML *pmlzPtr;
    const ModeHandler *sourcePtr;
    const RectangularGrid *gridPtr;


    void record_slice(const multi_array<std::complex<double>, 2> &buffer,
                      multi_array<std::complex<double>, 2> &storage, int idx, bool slice_y) const;



    /**
     * Do one step in propagation direction with the Crank-Nicolson scheme
     * @param field old field value. This value is not modified.
     * @param x position along the propagation direction
     * @param dx step size in propagation direction
     * @return new field value
     */
    [[nodiscard]] multi_array<std::complex<double>, 2> do_step_cn(const multi_array<std::complex<double>, 2> &field, double x,
                                                    double dx) const;
};


#endif //SOLVER_H
