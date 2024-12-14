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
#include <boost/multi_array.hpp>


/**
 * Geometry class. It can currently only do extrusion in the x direction, so a base angle of 90 degrees.
 */
class Geometry {
public:
    Geometry(const std::vector<Shape> &shapes, const std::unordered_map<std::string, XS> &xs_map);

    /**
     * Get index value in a single point (x,y,z). If you need the index in a line or a plane, it is recommended to use
     * the corresponding method for it as it is MUCH faster.
     * @return index value
     */
    [[nodiscard]] double get_index(double x, double y, double z) const;

    /**
     * Get index along a line in the z dimension. This saves time compared to looping over get_index(),
     * because it uses the same stack for all points in the array. If you need a slice, it is recommended to use the
     * even faster get_index_plane() method.
     * @param x x position
     * @param y y position
     * @param zgrid z positions
     * @return index values
     */
    std::vector<double> get_index_line(double x, double y, const std::vector<double>& zgrid) const;

    /**
     * Get index in a plane with normal direction the x direction. This maximally re-uses cached information.
     * So if you need a slice of the index, it is strongly recommended to use this method
     * instead of looping over get_index().
     * @param x x position
     * @param ygrid y positions
     * @param zgrid z positions
     * @return index values
     */
    boost::multi_array<double, 2> get_index_plane(double x, const std::vector<double> &ygrid, const std::vector<double>& zgrid) const;

private:
    std::vector<Shape> m_shapes;
    const std::unordered_map<std::string, XS> xs_map;

    XS get_xs(double x, double y) const;
};


#endif //GEOMETRY_H
