#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../src/TriDiag.h"
#include "../src/Polygon.h"
#include "../src/AuxiliaryFunctions.h"
#include "../src/Geometry.h"
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

TEST_CASE("Linspace test") {
    std::vector<double> my_range = AuxiliaryFunctions::linspace(-13.5, 26.0, 19);
    CHECK(my_range.front() == -13.5);
    CHECK(my_range[9] == 6.25);
    CHECK(my_range.back() == 26.0);
}

TEST_CASE("Check geometry class") {
    // we have a few overlapping material blocks here.
    constexpr double back_index = 1.4;
    constexpr double index1 = 2.0;
    constexpr double index2 = 1.7;
    constexpr double index3 = 2.4;
    const Polygon polygon1({Point(-1.0, 0.0), Point(0.3, 0.0), Point(0.3, 1.0), Point(0.0, 1.0)});
    const Polygon polygon2({Point(0.3, 0.0), Point(0.75, 0.0), Point(0.75, 0.65), Point(0.3, 0.65)});
    const Polygon polygon3({Point(-0.2, -0.4), Point(0.5, -0.4), Point(0.5, 0.4), Point(-0.2, 0.4)});

    const Shape shape1(polygon1, -0.5, 0.5, index1);
    const Shape shape2(polygon2, -0.4, 1.1, index2);
    const Shape shape3(polygon3, -1.0, 1.0, index3);
    const std::vector shapes = {shape1, shape2, shape3};
    Geometry geometry(shapes, back_index);

    CHECK(geometry.get_index(-2.0,0,0)==back_index); //outside of all shapes
    CHECK(geometry.get_index(0,0,0)==index3); //should be in the third shape
    CHECK(geometry.get_index(-0.5,0.4,0)==index1);
    CHECK(geometry.get_index(0.6,0.2,0)==index2);

    //overlapping here, but shape3 takes precedence because it was added later
    CHECK(geometry.get_index(0,0.16,0)==index3);
    CHECK(geometry.get_index(0.37,0.2,0)==index3);
    CHECK(geometry.get_index(0.37,0.2,1.4)==back_index);
    CHECK(geometry.get_index(0.37,0.2,1.05)==index2);
}
