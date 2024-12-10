//
// Created by mike on 12/9/24.
//

#include "Layer.h"

Layer::Layer(const double zmin, const double zmax, const double index) : index(index),
                                                                         zmin(zmin), zmax(zmax) {
}

bool Layer::is_in_layer(const double z) const {
    return z >= zmin && z <= zmax;
}

double Layer::get_index() const {
    return index;
}
