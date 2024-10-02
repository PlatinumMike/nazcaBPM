//
// Created by mike on 9/29/24.
//

#ifndef ENGINE_H
#define ENGINE_H

#include <complex>
#include <vector>

#include "Parameters.h"

/**
* Note, the coordinate system is the same as in Fig. 2.1 of the book "Beam propagation method for design of optical waveguide devices", G. Pedrola, 2016.
* So z is the propagation direction, x is normal out of the PIC surface, and y lies tangent to the PIC. This forms a right-handed coordinate system.
* It's a bit awkward because TE mode points in y direction, and TM in x. But anyway, let's be consistent with the book.
*/
class Engine {
public:
    explicit Engine(const std::string &fileName);

    /**
     * Run the actual simulation
     */
    void run();

    double get_refractive_index(double x, double y, double z) const;

    double get_min_dz(double dx, double dy) const;

    std::vector<std::complex<double> > get_initial_profile(const std::vector<double> &xgrid,
                                                           const std::vector<double> &ygrid) const;

    //derivative of the field vector w.r.t. z, so du/dz.
    std::vector<std::complex<double> > get_derivative(const std::vector<std::complex<double> > &field,
                                                      const std::vector<double> &xgrid,
                                                      const std::vector<double> &ygrid,
                                                      double z) const;

    // Update field with one step of the forward Euler scheme
    std::vector<std::complex<double> > do_step_euler(const std::vector<std::complex<double> > &field,
                                                     const std::vector<double> &xgrid,
                                                     const std::vector<double> &ygrid, double z, double dz) const;

    // Update field with one step of the explicit fourth order Runge-Kutta scheme.
    // Using Runge Kutta 4 for explicit time stepping. This has a limited dz step, proportional to dx^2.
    // Which may make it slower than the Crank-Nickelson scheme, it needs to do many more steps.
    // but each individual step is faster, also it has a smaller error per step than CN. So not clear which is faster, tbd.
    std::vector<std::complex<double> > do_step_rk4(const std::vector<std::complex<double> > &field,
                                                   const std::vector<double> &xgrid,
                                                   const std::vector<double> &ygrid,
                                                   double z, double dz) const;

    //get conductivity for the pml, in units of omega*eps0
    double get_conductivityx(double x) const;

    double get_conductivityy(double y) const;

    std::complex<double> get_qfactorx(double x, double y, double z) const;

    std::complex<double> get_qfactory(double x, double y, double z) const;

    std::vector<double> get_intensity(const std::vector<std::complex<double> > &field);

private:
    Parameters _inputs;
    // TODO these values do not need to be accessible globally in Engine, pass them along.
    int numx;
    int numy;
    int numz;
    double domain_len_x;
    double domain_len_y;
    double beta_ref;
    double k0;
    double max_index;
    double reference_index;
    double pml_thickness;
    double pml_strength;
    Geometry *geometry;


    // helper function to get the conductivity
    double get_conductivity_base(double x, double xmin, double xmax) const;

    // array indexing helper
    int get_neighbour_index(int i, int j, int ny, int direction = 0) const;

    void record_slice(const std::vector<std::complex<double> > &buffer,
                      std::vector<std::complex<double> > &storage) const;

    int dump_index_slice(const std::string &filename, char direction, double slice_position,
                         const std::vector<double> &xgrid,
                         const std::vector<double> &ygrid, const std::vector<double> &zgrid) const;
};


#endif //ENGINE_H
