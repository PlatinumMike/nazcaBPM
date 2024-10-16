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

std::complex<double> PML::get_pml_factor_derivative(const double x, const double pml_index) const {
    const auto factor = get_pml_factor(x, pml_index);
    return std::complex<double>{0.0, get_conductivity_derivative(x)} * factor * factor / (pml_index * pml_index);
}

double PML::get_conductivity(const double x) const {
    const double x_left = x_min + pml_thickness;
    const double x_right = x_max - pml_thickness;
    //second order polynomial PML strength
    if (x < x_left) {
        const double base_value = (x - x_left) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else if (x > x_right) {
        const double base_value = (x - x_right) / pml_thickness;
        return pml_strength * base_value * base_value;
    } else {
        return 0.0;
    }
}

double PML::get_conductivity_derivative(const double x) const {
    const double x_left = x_min + pml_thickness;
    const double x_right = x_max - pml_thickness;
    if (x < x_left) {
        const double base_value = (x - x_left) / pml_thickness;
        return pml_strength * base_value * 2.0 / pml_thickness;
    } else if (x > x_right) {
        const double base_value = (x - x_right) / pml_thickness;
        return pml_strength * base_value * 2.0 / pml_thickness;
    } else {
        return 0.0;
    }
}
