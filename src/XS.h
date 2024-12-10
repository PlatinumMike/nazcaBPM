//
// Created by mike on 12/9/24.
//

#ifndef XS_H
#define XS_H

#include <vector>
#include "Layer.h"

/**
 * Class to describe the cross-section (XS) of a given platform.
 * This lists all the
 */
class XS {
public:
    explicit XS(double background_index);

    void append_layer(const Layer &layer);

    [[nodiscard]] double get_index(double z, bool is_inside_core) const;

private:
    const double background_index;
    std::vector<Layer> core_layers;
    std::vector<Layer> substrate_layers;
};


#endif //XS_H
