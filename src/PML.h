//
// Created by mike on 10/11/24.
//

#ifndef PML_H
#define PML_H

#include <complex>

// class to implement a perfectly matched layer (PML) with a parabolic conductivity profile.

class PML {
public:
    PML(double thickness, double strength, double xmin, double xmax);

    /**
     * correction term to describe PML
     * @param x coordinate x
     * @param y coordinate y
     * @param z coordinate z
     * @param pml_index refractive index at the point (x,y,z). Just choose this equal to the index value of the underlying medium.
     * @return
     */
    [[nodiscard]] std::complex<double> get_pml_factor(double x, double y, double z, double pml_index) const;

private:
    double pml_thickness;
    double pml_strength;
    double x_min; //minimal coordinate value. Writen as x here, but it can represent any of the (x,y,z) coordinates.
    double x_max; //max coordinate value.

    //get conductivity for the pml, in units of omega*eps0
    [[nodiscard]] double get_conductivity(double x) const;

    // helper function to get the conductivity
    [[nodiscard]] double get_conductivity_base(double x, double xmin, double xmax) const;
};


#endif //PML_H
