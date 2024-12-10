//
// Created by mike on 10/1/24.
//

#include "Geometry.h"

#include <ranges>

Geometry::Geometry(const std::vector<Shape> &shapes,
                   const std::unordered_map<std::string, XS> &xs_map) : xs_map(xs_map) {
    m_shapes = shapes;
}

double Geometry::get_index(const double x, const double y, const double z) const {
    // loop over all shapes, search for a hit of x,y,z.
    for (const auto &shape: std::ranges::views::reverse(m_shapes)) {
        if (shape.point_in_shape(x, y)) {
            const auto xs = xs_map.at(shape.get_xs_name());
            return xs.get_index(z, true);
        }
    }

    //list exhausted, point is in cladding or substrate.
    return xs_map.at("default").get_index(z, false);
}
