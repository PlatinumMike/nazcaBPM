//
// Created by mike on 12/9/24.
//

#include "XS.h"

XS::XS(const double background_index) : background_index(background_index) {
}

void XS::append_layer(const Layer &layer) {
    if (layer.is_substrate()) {
        substrate_layers.push_back(layer);
    } else {
        core_layers.push_back(layer);
    }
}

double XS::get_index(const double z, const bool is_inside_core) const {
    //check if point is inside the substrate layers
    for (const auto &layer: substrate_layers) {
        if (layer.is_in_layer(z)) {
            return layer.get_index();
        }
    }
    // if the xy are inside the polygon that describes the core shape, also check those layers
    if (is_inside_core) {
        for (const auto &layer: core_layers) {
            if (layer.is_in_layer(z)) {
                return layer.get_index();
            }
        }
    }

    // The point is outside any of the specified layers, so return the background index
    return background_index;
}
