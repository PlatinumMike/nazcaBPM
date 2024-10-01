#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../src/TriDiag.h"
#include "../src/Polygon.h"
#include <cmath>

TEST_CASE("Tri-diagonal matrix solve") {
    //use a simple 4x4 matrix

    std::vector<double> diagonal = {1.0, -1.2, 0.5, 6.6};
    std::vector<double> upper = {3.33, 1.22, -9.0};
    std::vector<double> lower = {2.111, 1.0, 1.0};
    std::vector<double> solution_exact = {1.0, 0.0, -0.5, -0.13};
    std::vector<double> solution_thomas;


    auto size = solution_exact.size();
    std::vector<double> rhs(size);
    rhs[0] = diagonal[0] * solution_exact[0] + upper[0] * solution_exact[1];
    rhs[size - 1] = diagonal[size - 1] * solution_exact[size - 1] + lower[size - 2] * solution_exact[size - 2];
    for (int i = 1; i < size - 1; i++) {
        rhs[i] = lower[i - 1] * solution_exact[i - 1] + diagonal[i] * solution_exact[i] + upper[i] * solution_exact[
                     i + 1];
    }

    solution_thomas = TriDiag::solve_thomas(lower, diagonal, upper, rhs);

    double rms_error = 0.0;
    for (size_t i = 0; i < size; i++) {
        double error = solution_thomas[i] - solution_exact[i];
        rms_error += error * error;
    }
    rms_error /= static_cast<double>(solution_thomas.size());
    rms_error = std::sqrt(rms_error);


    // it will never be exact due to rounding errors, but if it is within tolerance that is ok.
    double tolerance = 1.0e-14;
    CHECK(rms_error < tolerance);
}

TEST_CASE("Check ray-casting polygon") {
    //example shape in klayout
    std::vector<Point> poly_points = {
        {32.61000, -25.91100},
        {29.04400, -21.92500},
        {32.40100, -17.39900},
        {21.16100, -20.45600},
        {16.81500, -17.63900},
        {31.44100, -11.58400},
        {45.13900, -19.58700}
    };
    Polygon polygon(poly_points);

    Point test_point1 = {0, 0}; //out
    Point test_point2 = {31, -18}; //out
    Point test_point3 = {38, -20}; //in
    Point test_point4 = {25.5, -15}; //in

    CHECK(!polygon.point_inside_polygon(test_point1));
    CHECK(!polygon.point_inside_polygon(test_point2));
    CHECK(polygon.point_inside_polygon(test_point3));
    CHECK(polygon.point_inside_polygon(test_point4));
}
