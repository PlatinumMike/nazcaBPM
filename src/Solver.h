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
using boost::multi_array;

// todo: for the current implementation many arrays are copied at every iteration.
// Performance can be improved by overwriting them in place, at least where that is safe.


class Solver {
public:
    Solver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source, double xmin,
           double xmax, double ymin, double ymax, double zmin, double zmax, int numx, int numy, int numz,
           double scheme_parameter, double k0, double reference_index);

    void run();

    multi_array<double, 2> get_intensity(const multi_array<std::complex<double>, 2> &field) const;

private:
    //default values grid
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    double zmin = 0.0;
    double zmax = 0.0;
    int numx = 0;
    int numy = 0;
    int numz = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double scheme_parameter = 0.5;
    double k0;
    double reference_index;
    // pointers to simulation objects
    const Geometry *geometryPtr;
    const PML *pmlxPtr;
    const PML *pmlyPtr;
    const ModeHandler *sourcePtr;

    [[nodiscard]] double get_refractive_index(double x, double y, double z) const;

    void record_slice(const multi_array<std::complex<double>, 2> &buffer,
                      multi_array<std::complex<double>, 2> &storage, int idz, bool slice_x) const;

    void dump_index_slice(const std::string &filename, char direction, double slice_position,
                          const std::vector<double> &grid_coordinate1,
                          const std::vector<double> &grid_coordinate2) const;

    multi_array<double, 1> vector_to_multi_array(const std::vector<double> &vec) const;

    /**
     * Do one step in propagation direction with the Crank-Nicolson scheme
     * @param field old field value. This value is not modified.
     * @param xgrid grid in x direction
     * @param ygrid grid in y direction
     * @param z position along the propagation direction
     * @param dz step size in propagation direction
     * @return new field value
     */
    multi_array<std::complex<double>, 2> do_step_cn(const multi_array<std::complex<double>, 2> &field,
                                                    const std::vector<double> &xgrid, const std::vector<double> &ygrid,
                                                    double z,
                                                    double dz) const;

    multi_array<std::complex<double>, 2> get_rhs(const multi_array<std::complex<double>, 2> &field,
                                                 const std::vector<double> &xgrid,
                                                 const std::vector<double> &ygrid, double z, double dz) const;

    //applies the operator (1+prefactorGx) or (1+prefactorGy). You can decide which one to do based on the neighbor values you pass into the function.
    std::complex<double> apply_right_hand_operator(double xmid, double xneighbor_previous, double xneighbor_next,
                                                   double ymid,
                                                   double yneighbor_previous, double yneighbor_next,
                                                   double neighbor_distance,
                                                   double z,
                                                   std::complex<double> preFactor,
                                                   std::complex<double> fieldValueMid,
                                                   std::complex<double> fieldValue_previous,
                                                   std::complex<double> fieldValue_next, const PML *pmlPtr) const;

    std::vector<std::complex<double> > solve_system(const std::vector<double> &xmid,
                                                    const std::vector<double> &xneighbor_previous,
                                                    const std::vector<double> &xneighbor_next,
                                                    const std::vector<double> &ymid,
                                                    const std::vector<double> &yneighbor_previous,
                                                    const std::vector<double> &yneighbor_next,
                                                    double neighbor_distance,
                                                    double z, std::complex<double> preFactor,
                                                    const std::vector<std::complex<double> > &rhs_slice,
                                                    const PML *pmlPtr) const;
};


#endif //SOLVER_H
