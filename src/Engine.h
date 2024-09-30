//
// Created by mike on 9/29/24.
//

#ifndef ENGINE_H
#define ENGINE_H

#include <complex>
#include <numbers>
#include <vector>

/**
* Note, the coordinate system is the same as in Fig. 2.1 of the book "Beam propagation method for design of optical waveguide devices", G. Pedrola, 2016.
* So z is the propagation direction, x is normal out of the PIC surface, and y lies tangent to the PIC. This forms a right-handed coordinate system.
* It's a bit awkward because TE mode points in y direction, and TM in x. But anyway, let's be consistent with the book.
*/
class Engine {
public:
    Engine();

    /**
     * Run the actual simulation
     */
    void run();

    double get_refractive_index(double x, double y, double z) const;

    double get_min_dz(double dx,  double dy) const;

    std::vector<std::complex<double> > get_initial_profile(const double *xgrid, const double *ygrid) const;

    //derivative of the field vector w.r.t. z, so du/dz.
    std::vector<std::complex<double> > get_derivative(const std::vector<std::complex<double> > &field, double *xgrid,
                                                      double *ygrid,
                                                      double z) const;

    // Update field with one step of the forward Euler scheme
    std::vector<std::complex<double> > do_step_euler(const std::vector<std::complex<double> > &field, double *xgrid,
                                                     double *ygrid, double z, double dz) const;

    // Update field with one step of the explicit fourth order Runge-Kutta scheme.
    // Using Runge Kutta 4 for explicit time stepping. This has a limited dz step, proportional to dx^2.
    // Which may make it slower than the Crank-Nickelson scheme, it needs to do many more steps.
    // but each individual step is faster, also it has a smaller error per step than CN. So not clear which is faster, tbd.
    std::vector<std::complex<double> > do_step_rk4(const std::vector<std::complex<double> > &field, double *xgrid,
                                                   double *ygrid,
                                                   double z, double dz) const;

    //get conductivity for the pml, in units of omega*eps0
    double get_conductivityx(double x) const;

    double get_conductivityy(double y) const;

    std::complex<double> get_qfactorx(double x, double y, double z) const;

    std::complex<double> get_qfactory(double x, double y, double z) const;

    std::vector<double> get_intensity(const std::vector<std::complex<double> >& field);

private:
    //TODO: make these inputs that are read from a file
    double background_index = 1.5; //background refractive index
    double core_index = 2.0;
    double wl = 1.55; //vacuum wavelength (um)
    int resolution_x = 10; //number of gridpoints per unit length in x direction
    int resolution_y = 10; //number of gridpoints per unit length in y direction
    int numx = 1; //number of gridpoints
    int numy = 1;
    int numz = 1;
    double domain_len_x = 4; //length of the domain in x direction
    double domain_len_y = 6.0; //length of the domain in y direction
    double domain_len_z = 10.0; //length of the domain in z direction
    double reference_index = 1.6;
    double k0 = 2 * std::numbers::pi / wl; //vacuum wavenumber
    double beta_ref = k0 * reference_index;

    double pml_strength = 5.0;
    double pml_thickness = 1.0; //thickness of the PML (um)

    // helper function to get the conductivity
    double get_conductivity_base(double x, double xmin, double xmax) const;

    // array indexing helper
    int get_neighbour_index(int i, int j, int ny, int direction = 0) const;

    void record_slice(const std::vector<std::complex<double> > &buffer,std::vector<std::complex<double> > &storage);

};


#endif //ENGINE_H
