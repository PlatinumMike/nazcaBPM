//
// Created by mike on 9/29/24.
//

#include "AuxiliaryFunctions.h"

void AuxiliaryFunctions::linspace(double *array, const double start, const double stop, const int count){
    double delta = (stop - start) / (count - 1);
    for (int i = 0; i < count; ++i) {
        array[i] = start + delta * i;
    }
}
