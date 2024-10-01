//
// Created by mike on 10/1/24.
//

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include "Polygon.h"


// shapes in nazca are described by their polygon, then extruded vertically. Also we will assume they have a uniform refractive index.
// You could
struct Shape {
    Polygon poly;
    double zmin = 0;
    double zmax = 0;
    double refractive_index = 1.0;
};


/**
 * Geometry class. It can currently only do extrusion in the x direction, so a base angle of 90 degrees.
 */
class Geometry {
public:
    Geometry(std::vector<Shape> shapes, double background_index);

    //loop over all shapes, find if the point is inside, and return the index.
    double get_index(double x, double y, double z);

private:
    // todo: keep track of a temporary sub-list for each z slice for only the shapes that are present in that slice.
    // this speeds up the searching.
    std::vector<Shape> m_shapes;
    double m_background_index;

    bool point_in_shape(double x, double y, double z, const Shape& shape);
};


#endif //GEOMETRY_H
