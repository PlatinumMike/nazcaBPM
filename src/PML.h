//
// Created by mike on 10/11/24.
//

#ifndef PML_H
#define PML_H

#include <complex>

// class to implement a perfectly matched layer (PML) with a parabolic conductivity profile.

class PML {
public:
    /**
     * creates a PML object that defines a PML near xmin and xmax.
     * E.g. if the x domain spans [0,10], and the thickness is 1, a PML will be applied from 0 to 1, and 9 to 10.
     * The conductivity profile is zero at the entrance of the PML and ramps up parabolically towards the domain boundaries.
     * @param thickness PML thickness in length units
     * @param strength PML strength in units of omega*eps0
     * @param xmin min coordinate value of the domain.
     * @param xmax max coordinate value of the domain.
     */
    PML(double thickness, double strength, double xmin, double xmax);

    /**
     * correction term to describe PML
     * @param x coordinate normal to the PML.
     * @param pml_index refractive index at the point of interest. Just choose this equal to the index value of the underlying medium.
     * @return pml factor
     */
    [[nodiscard]] std::complex<double> get_pml_factor(double x, double pml_index) const;

    /**
     * Derivative w.r.t. x of correction term to describe PML
     * not to be used in the simulation, this is purely for testing purposes.
     * @param x coordinate normal to the PML.
     * @param pml_index refractive index at the point of interest. Just choose this equal to the index value of the underlying medium.
     * @return derivative of pml factor
     */
    [[nodiscard]] std::complex<double> get_pml_factor_derivative(double x, double pml_index) const;

    //get conductivity for the pml, in units of omega*eps0
    [[nodiscard]] double get_conductivity(double x) const;

private:
    double pml_thickness;
    double pml_strength;
    double x_min; //minimal coordinate value. Writen as x here, but it can represent any of the (x,y,z) coordinates.
    double x_max; //max coordinate value.

    //get dervivative of conductivity for the pml, w.r.t. x.
    [[nodiscard]] double get_conductivity_derivative(double x) const;
};

#endif //PML_H
