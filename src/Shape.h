//
// Created by mike on 12/9/24.
//

#ifndef SHAPE_H
#define SHAPE_H
#include <string>

#include "Polygon.h"

// shapes in nazca are described by their polygon, then extruded vertically. Also we will assume they have a uniform refractive index.
// you could embed the extra z and index information into the polygon, but this is kept apart on purpose. This way the nazca information is kept separate.
// the ray casting only needs to be done for the x,y coordinates, it does not need any other information.


class Shape : public Polygon {
public:
    Shape(const std::vector<Point> &points_input, std::string xs_name);

    /**
     * Check if a point (x,y) is inside the polygon of the Shape
     * @param x x coordinate
     * @param y y coordinate
     * @return true if (x,y) is in the shape
     */
    [[nodiscard]] bool point_in_shape(double x, double y) const;

    [[nodiscard]] std::string get_xs_name() const;

private:
    // each shape can have it's own cross section.
    // instead of storing a 100 copies of the xs, we just store the key to lookup the xs later in a map of cross sections.
    std::string xs_name;
};


#endif //SHAPE_H
