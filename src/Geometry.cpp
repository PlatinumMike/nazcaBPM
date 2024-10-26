//
// Created by mike on 10/1/24.
//

#include "Geometry.h"

#include <ranges>

Geometry::Geometry(const std::vector<Shape> &shapes, const double background_index) : m_background_index(background_index) {
    m_shapes = shapes;
}

double Geometry::get_index(const double x, const double y, const double z) const {
    // loop over all shapes, search for a hit of x,y,z.
    for (const auto &shape: std::ranges::views::reverse(m_shapes)) {
        if (point_in_shape(x, y, z, shape)) {
            return shape.refractive_index;
        }
    }

    //list exhausted, return background index.
    return m_background_index;
}

bool Geometry::point_in_shape(const double x, const double y, const double z, const Shape &shape) const {
    if (z < shape.zmin || z > shape.zmax) {
        //definitely not inside the shape.
        return false;
    } else {
        //maybe in the shape, search polygon
        const Point point = {x, y};
        return shape.poly.point_inside_polygon(point);
    }
}
