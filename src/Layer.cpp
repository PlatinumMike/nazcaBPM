//
// Created by mike on 12/9/24.
//

#include "Layer.h"

Layer::Layer(const double zmin, const double zmax, const double index, const bool isSubstrate) : index(index),
    zmin(zmin), zmax(zmax), isSubstrate(isSubstrate) {
}

bool Layer::is_in_layer(const double z) const {
    return z >= zmin && z <= zmax;
}

double Layer::get_index() const {
    return index;
}

bool Layer::is_substrate() const {
    return isSubstrate;
}
