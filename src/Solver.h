//
// Created by mike on 10/11/24.
//

#ifndef SOLVER_H
#define SOLVER_H

#include <complex>
#include "Geometry.h"
#include <boost/multi_array.hpp>


#include "PML.h"
#include "RectangularGrid.h"
using boost::multi_array;

// todo: for the current implementation many arrays are copied at every iteration.
// Performance can be improved by overwriting them in place, at least where that is safe.


class Solver {
public:
    Solver(const Geometry &geometry, const PML &pmly, const PML &pmlz,
           const RectangularGrid &grid,
           double scheme_parameter, double k0, double reference_index);


    multi_array<double, 2> get_intensity(const multi_array<std::complex<double>, 2> &field) const;

    /**
     * Interpolate the internal field on a given new grid.
     * @param ygrid_new new y coordinates
     * @param zgrid_new new z coordinates
     * @return field on the new coordinates.
     */
    multi_array<std::complex<double>, 2> interpolate_field(const std::vector<double> &ygrid_new,
                                                           const std::vector<double> &zgrid_new) const;

protected:
    /**
     * Do one step in propagation direction with the Crank-Nicolson scheme
     * @param field old field value. This value is not modified.
     * @param x position along the propagation direction
     * @param dx step size in propagation direction
     * @param propagation_factor i*dx / (2 * k0 * reference_index);
     * @return new field value
     */
    [[nodiscard]] multi_array<std::complex<double>, 2> do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                                  double x,
                                                                  double dx,
                                                                  std::complex<double> propagation_factor) const;

    const RectangularGrid *gridPtr;
    multi_array<std::complex<double>, 2> internal_field;
    const double k0;
    const double reference_index;

private:
    const double scheme_parameter = 0.5;
    // pointers to simulation objects
    const Geometry *geometryPtr;
    const PML *pmlyPtr;
    const PML *pmlzPtr;
};


#endif //SOLVER_H
