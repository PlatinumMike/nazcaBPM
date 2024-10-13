//
// Created by mike on 10/11/24.
//

#include "PML.h"

PML::PML(const double thickness, const double strength, const double xmin,
         const double xmax) : pml_thickness(thickness),
                              pml_strength(strength), x_min(xmin),
                              x_max(xmax) {
}

std::complex<double> PML::get_pml_factor(const double x, const double pml_index) const {
    return 1.0 / std::complex{1.0, -get_conductivity(x) / (pml_index * pml_index)};
}

double PML::get_conductivity(const double x) const {
    //second order polynomial PML strength
    if (x < x_min + pml_thickness) {
        const double base_value = ((x_min + pml_thickness) - x) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else if (x > x_max - pml_thickness) {
        const double base_value = (x - (x_max - pml_thickness)) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else {
        return 0.0;
    }
}
