//
// Created by mike on 10/1/24.
//

#ifndef POLYGON_H
#define POLYGON_H
#include <vector>

/**
 * Polygon class. This uses nazca/klayout style (x,y) coordinates. The z direction is the normal direction of the chip.
 * So out of the plane. These coordinates are not the same as the BPM coordinates.
 * x_nazca = z_bpm
 * y_nazca = -y_bpm
 * z_nazca = x_bpm
 * Polygons in the polygon class can be described with a vector of points (x_nazca, y_nazca).
 */
struct Point {
    double x, y;
};

class Polygon {
public:
    // Note, for a polygon with N vertices, you have to pass in a vector of size N.
    // in some conventions people add the starting point again at the end, to "close" the polygon.
    // Do not do this! This method expects unique points only!
    explicit Polygon(std::vector<Point> points_input);

    // evaluate if a point is inside the polygon.
    [[nodiscard]] bool point_inside_polygon(Point point) const;

private:
    std::vector<Point> points;
};


#endif //POLYGON_H
