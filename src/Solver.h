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


class Solver {
public:
    Solver(const Geometry &geometry, const PML &pmlx, const PML &pmly, const ModeHandler &source, double xmin,
           double xmax, double ymin, double ymax, double zmin, double zmax, int numx, int numy, int numz,
           double scheme_parameter);

    void run();

    //derivative of the field vector w.r.t. z, so du/dz.
    multi_array<std::complex<double>, 2> get_derivative(const multi_array<std::complex<double>, 2> &field,
                                                        const std::vector<double> &xgrid,
                                                        const std::vector<double> &ygrid,
                                                        double z) const;


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
    // pointers to simulation objects
    const Geometry *geometryPtr;
    const PML *pmlxPtr;
    const PML *pmlyPtr;
    const ModeHandler *sourcePtr;

    double get_refractive_index(double x, double y, double z) const;

    void record_slice(const multi_array<std::complex<double>, 2> &buffer,
                      multi_array<std::complex<double>, 2> &storage, int idz) const;

    void dump_index_slice(const std::string &filename, char direction, double slice_position,
                          const std::vector<double> &grid_coordinate1,
                          const std::vector<double> &grid_coordinate2) const;

    multi_array<double, 1> vector_to_multi_array(const std::vector<double> &vec) const;

    // one step with the Crank-Nicolson scheme
    multi_array<std::complex<double>, 2> do_step_cn(multi_array<std::complex<double>, 2> field,
                                                    std::vector<double> xgrid, std::vector<double> ygrid, double z,
                                                    double dz);
};


#endif //SOLVER_H
