//
// Created by mike on 10/1/24.
//

#include "Polygon.h"

#include <utility>

Polygon::Polygon(std::vector<Point> points_input) {
    points = std::move(points_input);
}

//source: adapted from: https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
bool Polygon::point_inside_polygon(const Point point) const {
    const int num_vertices = static_cast<int>(points.size());
    const double x = point.x;
    const double y = point.y;
    bool inside = false;

    // Store the first point in the polygon and initialize
    // the second point
    Point p1 = points[0], p2{};

    // Loop through each edge in the polygon
    for (int i = 1; i <= num_vertices; i++) {
        // Get the next point in the polygon
        p2 = points[i % num_vertices];

        // Check if the point is above the minimum y
        // coordinate of the edge
        if (y > std::min(p1.y, p2.y)) {
            // Check if the point is below the maximum y
            // coordinate of the edge
            if (y <= std::max(p1.y, p2.y)) {
                // Check if the point is to the left of the
                // maximum x coordinate of the edge
                if (x <= std::max(p1.x, p2.x)) {
                    // Calculate the x-intersection of the
                    // line connecting the point to the edge
                    const double x_intersection
                            = (y - p1.y) * (p2.x - p1.x)
                              / (p2.y - p1.y)
                              + p1.x;

                    // Check if the point is on the same
                    // line as the edge or to the left of
                    // the x-intersection
                    if (p1.x == p2.x
                        || x <= x_intersection) {
                        // Flip the inside flag
                        inside = !inside;
                    }
                }
            }
        }

        // Store the current point as the first point for
        // the next iteration
        p1 = p2;
    }

    // Return the value of the inside flag
    return inside;
}
