//
// Created by mike on 12/9/24.
//

#ifndef LAYER_H
#define LAYER_H

/**
 * Class to describe cross-section layer information
 * Layers have a minimal and maximal z value, and a refractive index.
 * The refractive index is assumed to be uniform within the layer.
 */
class Layer {
public:
    Layer(double zmin, double zmax, double index);

    [[nodiscard]] bool is_in_layer(double z) const;

    [[nodiscard]] double get_index() const;

private:
    const double index;
    const double zmin;
    const double zmax;
};


#endif //LAYER_H
