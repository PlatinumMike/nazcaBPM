//
// Created by mike on 10/1/24.
//

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <vector>
#include <string>
#include <unordered_map>

#include "Shape.h"
#include "XS.h"


/**
 * Geometry class. It can currently only do extrusion in the x direction, so a base angle of 90 degrees.
 */
class Geometry {
public:
    Geometry(const std::vector<Shape> &shapes, const std::unordered_map<std::string, XS> &xs_map);

    //loop over all shapes, find if the point is inside, and return the index.
    [[nodiscard]] double get_index(double x, double y, double z) const;

private:
    // todo: keep track of a temporary sub-list for each x slice for only the shapes that are present in that slice.
    // this speeds up the searching.
    std::vector<Shape> m_shapes;
    const std::unordered_map<std::string, XS> xs_map;


};


#endif //GEOMETRY_H
