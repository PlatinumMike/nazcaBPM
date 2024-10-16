#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../src/TriDiag.h"
#include "../src/Polygon.h"
#include "../src/AuxiliaryFunctions.h"
#include "../src/Geometry.h"
#include "../src/GridInterpolator.h"
#include "../src/PML.h"
#include "../src/OperatorSuite.h"
#include <cmath>
#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;


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

    solution_thomas = TriDiag<double>::solve_thomas(lower, diagonal, upper, rhs);

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

TEST_CASE("Check grid interpolation") {
    int numx = 6;
    int numy = 4;
    auto xgrid = AuxiliaryFunctions::linspace(-4.0, 4.0, numx);
    auto ygrid = AuxiliaryFunctions::linspace(-2.0, 1.5, numy);
    boost::multi_array<double, 2> dataset(boost::extents[numx][numy]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            // data = 3x^2-0.3y+8
            dataset[i][j] = 3.0 * xgrid[i] * xgrid[i] - 0.3 * ygrid[j] + 8.0;
        }
    }
    double external_value = 0.0;


    // get value at some sample points
    GridInterpolator<double> grid_interpolator(xgrid, ygrid, dataset, 0.0);

    // outside values
    CHECK(grid_interpolator.get_value(0.1, 16.0)==external_value);
    CHECK(grid_interpolator.get_value(-5.0, 0.0)==external_value);

    // add some internal points in the middle
    CHECK(grid_interpolator.get_value(-1.12, 0.3)==doctest::Approx(12.902));
    CHECK(grid_interpolator.get_value(-1.82,0.45)==doctest::Approx(19.577));
    CHECK(grid_interpolator.get_value(2.0,0.11)==doctest::Approx(21.407));
    CHECK(grid_interpolator.get_value(2.2,0.11)==doctest::Approx(23.327));

    // corner values
    CHECK(grid_interpolator.get_value(4.0, 1.5)==doctest::Approx(55.55));
    CHECK(grid_interpolator.get_value(4.0, -2.0)==doctest::Approx(56.6));
    CHECK(grid_interpolator.get_value(-4.0, -2.0)==doctest::Approx(56.6));
    CHECK(grid_interpolator.get_value(-4.0, 1.5)==doctest::Approx(55.55));

    // edge values
    CHECK(grid_interpolator.get_value(-4.0, 0.0)==doctest::Approx(56.0));
    CHECK(grid_interpolator.get_value(4.0, 0)==doctest::Approx(56.0));

    // internal point which happens to be exactly on a grid point.
    CHECK(grid_interpolator.get_value(-2.4, 1.0/3.0)==doctest::Approx(25.18));
}

TEST_CASE("Check PML") {
    PML pml(1.5, 4.0, 0.0, 10.0);
    double index = 2.0;

    CHECK(pml.get_conductivity(0.0)==4.0);
    CHECK(pml.get_conductivity(10.0)==4.0);
    CHECK(pml.get_conductivity(1.5)==0.0);
    CHECK(pml.get_conductivity(8.5)==0.0);
    CHECK(pml.get_conductivity(9.0)==4.0/9.0);

    CHECK(pml.get_pml_factor(5.0,index)==std::complex<double>{1.0, 0.0});
    CHECK(pml.get_pml_factor(9.0,index)==std::complex<double>{0.9878048780487805, 0.1097560975609756});
    CHECK(pml.get_pml_factor(10.0,index)==std::complex<double>{0.5, 0.5});
}


TEST_CASE("Check Gx operator") {
    int numx = 1000;
    int numy = 120;
    std::complex<double> amplitude = {0.3, 0.5}; //u0
    std::complex<double> preFactor = {-0.1, 0.8}; //p
    auto xgrid = AuxiliaryFunctions::linspace(-4.0, 4.0, numx);
    auto ygrid = AuxiliaryFunctions::linspace(-2.0, 6.0, numy);
    double k0 = 1.0;
    double reference_index = 1.6;
    double index = 1.5; //using a uniform medium.
    PML pmlx(1.0, 5.0, xgrid.front(), xgrid.back());

    //Using some example field: u = u0*sin(k0*x)*sin(2*k0*y)
    multi_array<std::complex<double>, 2> field(extents[numx][numy]);
    multi_array<std::complex<double>, 2> field_derivative(extents[numx][numy]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            field[i][j] = amplitude * sin(k0 * xgrid[i]) * sin(2.0 * k0 * ygrid[j]);
            field_derivative[i][j] = amplitude * k0 * cos(k0 * xgrid[i]) * sin(2.0 * k0 * ygrid[j]);
        }
    }

    // (1+p*Gx)u = (1+p/2*k0^2*(n^2-n0^2))*u+p*eta*eta'*u'+p*eta^2*u''
    // where prime is the derivative w.r.t. x, and eta is the PML factor.
    multi_array<std::complex<double>, 2> reference_values(extents[numx][numy]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                reference_values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                auto eta = pmlx.get_pml_factor(xgrid[i], index);
                auto eta_prime = pmlx.get_pml_factor_derivative(xgrid[i], index);
                auto coef1 = 1.0 + 0.5 * preFactor * k0 * k0 * (index * index - reference_index * reference_index) -
                             preFactor * k0 * k0 * eta * eta;
                auto coef2 = preFactor * eta * eta_prime;
                reference_values[i][j] = coef1 * field[i][j] + coef2 * field_derivative[i][j];
            }
        }
    }

    //now apply the discretized form of the operator.
    multi_array<std::complex<double>, 2> values(extents[numx][numy]);
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            if (i == 0 || j == 0 || i == numx - 1 || j == numy - 1) {
                //ignore boundary values.
                values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // y value on the mid-point and neighbors is the same, as we only apply the derivative in x direction here.
                values[i][j] = OperatorSuite::apply_right_hand_operator(xgrid[i], xgrid[i - 1], xgrid[i + 1], index,
                                                                        index, index, reference_index, k0, preFactor,
                                                                        field[i][j], field[i - 1][j], field[i + 1][j],
                                                                        &pmlx);
            }
        }
    }

    double rms_error = 0.0;
    for (int i = 0; i < numx; i++) {
        for (int j = 0; j < numy; j++) {
            double base = std::abs(values[i][j] - reference_values[i][j]);
            rms_error += base * base;
        }
    }
    rms_error = std::sqrt(rms_error / (numx * numy));

    //it will never be perfectly equal anyway, because it's an approximation on a grid. But the error should be small.
    double abs_tol = 1.0e-3;
    CHECK(rms_error <abs_tol);
}
