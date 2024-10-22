//
// Created by mike on 10/1/24.
//

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include "Polygon.h"


// shapes in nazca are described by their polygon, then extruded vertically. Also we will assume they have a uniform refractive index.
// you could embed the extra z and index information into the polygon, but this is kept apart on purpose. This way the nazca information is kept separate.
// the ray casting only needs to be done for the x,y coordinates, it does not need any other information.
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
    Geometry(const std::vector<Shape> &shapes, double background_index);

    //loop over all shapes, find if the point is inside, and return the index.
    double get_index(double x, double y, double z) const;

private:
    // todo: keep track of a temporary sub-list for each z slice for only the shapes that are present in that slice.
    // this speeds up the searching.
    std::vector<Shape> m_shapes;
    double m_background_index;

    //this uses nazca coordinates
    bool point_in_shape(double x, double y, double z, const Shape &shape) const;
};


#endif //GEOMETRY_H
