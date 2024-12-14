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
    const XS xs = get_xs(x, y);
    return xs.get_index(z);
}

std::vector<double> Geometry::get_index_line(const double x, const double y, const std::vector<double> &zgrid) const {
    // first determine the XS, which is shared between all z points. So we do not need to check for other XS again.
    const XS xs = get_xs(x, y);
    std::vector<double> index(zgrid.size());
    for (size_t i = 0; i < zgrid.size(); ++i) {
        index[i] = xs.get_index(zgrid[i]);
    }
    return index;
}

boost::multi_array<double, 2> Geometry::get_index_plane(const double x, const std::vector<double> &ygrid,
                                                        const std::vector<double> &zgrid) const {
    //need to re-check the XS for every x,y point, so no speed-up vs get_index_line.
    int numy = static_cast<int>(ygrid.size());
    int numz = static_cast<int>(zgrid.size());

    // get index in slice
    boost::multi_array<double, 2> index(boost::extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        std::vector<double> index_line = get_index_line(x, ygrid[i], zgrid);
        for (int j = 0; j < numz; j++) {
            index[i][j] = index_line[j];
        }
    }
    return index;
}

XS Geometry::get_xs(const double x, const double y) const {
    // loop over all shapes, search for a hit of x,y
    //todo: this is a rather expensive operation, maybe we can cache the XS names on the x,y grid for more rapid access later on.
    for (const auto &shape: std::ranges::views::reverse(m_shapes)) {
        if (shape.point_in_shape(x, y)) {
            return xs_map.at(shape.get_xs_name());
        }
    }

    //list exhausted, point is in cladding or substrate.
    return xs_map.at("default");
}
