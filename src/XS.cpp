//
// Created by mike on 12/9/24.
//

#include "XS.h"

XS::XS(const double background_index) : background_index(background_index) {
}

void XS::append_layer(const Layer &layer) {
    layers.push_back(layer);
}

double XS::get_index(const double z) const {
    //check if point is inside the layers
    for (const auto &layer: layers) {
        if (layer.is_in_layer(z)) {
            return layer.get_index();
        }
    }
    // The point is outside any of the specified layers, so return the background index
    return background_index;
}
