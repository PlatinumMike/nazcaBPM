//
// Created by mike on 12/9/24.
//

#ifndef XS_H
#define XS_H

#include <vector>
#include "Layer.h"

/**
 * Class to describe the cross-section (XS) of a given platform.
 * This lists all the different layers that exist in it.
 * Note, this is purely a step-wise 1D profile of the refractive index.
 */
class XS {
public:
    explicit XS(double background_index);

    void append_layer(const Layer &layer);

    [[nodiscard]] double get_index(double z) const;

private:
    const double background_index;
    std::vector<Layer> layers;
};


#endif //XS_H
