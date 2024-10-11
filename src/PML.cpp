//
// Created by mike on 10/11/24.
//

#include "PML.h"

PML::PML(double thickness, double strength, double xmin, double xmax) : pml_thickness(thickness),
                                                                        pml_strength(strength), x_min(xmin),
                                                                        x_max(xmax) {
}

std::complex<double> PML::get_pml_factor(double x, double y, double z, double pml_index) const {
    return 1.0 / std::complex{1.0, -get_conductivity(x) / (pml_index * pml_index)};
}

double PML::get_conductivity(const double x) const {
    return get_conductivity_base(x, x_min, x_max);
}

double PML::get_conductivity_base(const double x, double xmin, double xmax) const {
    //second order polynomial PML strength
    if (x < xmin + pml_thickness) {
        const double base_value = ((xmin + pml_thickness) - x) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else if (x > xmax - pml_thickness) {
        const double base_value = (x - (xmax - pml_thickness)) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else {
        return 0.0;
    }
}
