//
// Created by mike on 12/9/24.
//

#include "Shape.h"

#include <utility>

Shape::Shape(const std::vector<Point> &points_input, std::string xs_name): Polygon(points_input),
                                                                           xs_name(std::move(xs_name)) {
}

bool Shape::point_in_shape(const double x, const double y) const {
    const Point point = {x, y};
    return this->point_inside_polygon(point);
}

std::string Shape::get_xs_name() const {
    return xs_name;
}
