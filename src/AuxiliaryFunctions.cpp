//
// Created by mike on 9/29/24.
//

#include "AuxiliaryFunctions.h"

std::vector<double> AuxiliaryFunctions::linspace(const double start, const double stop, const int count) {
    std::vector<double> array(count);
    const double delta = (stop - start) / (count - 1.0);
    for (int i = 0; i < count; ++i) {
        array[i] = start + delta * i;
    }
    return array;
}
