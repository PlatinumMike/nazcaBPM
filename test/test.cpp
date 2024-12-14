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

TEST_CASE("Layers") {
    Layer layer(-0.1, 2.8, 1.99);
    CHECK(layer.get_index() == 1.99);
    CHECK(layer.is_in_layer(0.0));
    CHECK(layer.is_in_layer(2.8));
    CHECK(!layer.is_in_layer(-1.0));
}

TEST_CASE("Cross-section") {
    constexpr double background_index = 1.45;
    XS xs(background_index);
    // no layers added yet, so it has to return the background index
    CHECK(xs.get_index(0.0)==background_index);
    CHECK(xs.get_index(1.0)==background_index);
    Layer layer1(0.0, 2.8, 1.99);
    Layer layer2(2.8, 4.0, 1.85);
    Layer layer3(4.0, 5.6, 1.22);
    xs.append_layer(layer1);
    xs.append_layer(layer2);
    xs.append_layer(layer3);
    CHECK(xs.get_index(-1.0)==background_index);
    CHECK(xs.get_index(0.0)==layer1.get_index());
    CHECK(xs.get_index(1.0)==layer1.get_index());
    CHECK(xs.get_index(2.8)==layer1.get_index());
    CHECK(xs.get_index(3.2)==layer2.get_index());
    CHECK(xs.get_index(4.1)==layer3.get_index());
    CHECK(xs.get_index(10.0)==background_index);
}

TEST_CASE("Check geometry class") {
    // we have a few overlapping material blocks here.
    constexpr double back_index = 1.4;
    constexpr double index1 = 2.0;
    constexpr double index2 = 1.7;
    constexpr double index3 = 2.4;
    // see "geometry_test.gds" for these shapes
    const std::vector<Point> points1 = {{-1.0, 0.0}, {0.3, 0.0}, {0.3, 1.0}, {0.0, 1.0}};
    const std::vector<Point> points2 = {{0.3, 0.0}, {0.75, 0.0}, {0.75, 0.65}, {0.3, 0.65}};
    const std::vector<Point> points3 = {{-0.2, -0.4}, {0.5, -0.4}, {0.5, 0.4}, {-0.2, 0.4}};

    // we give each polygon it's own distinct XS with just one layer
    Layer layer1(-0.5, 0.5, index1);
    Layer layer2(-0.4, 1.1, index2);
    Layer layer3(-1.0, 1.0, index3);
    XS xs1(back_index);
    xs1.append_layer(layer1);
    XS xs2(back_index);
    xs2.append_layer(layer2);
    XS xs3(back_index);
    xs3.append_layer(layer3);
    XS xs_default(back_index);
    std::unordered_map<std::string, XS> xs_map = {{"xs1", xs1}, {"xs2", xs2}, {"xs3", xs3}, {"default", xs_default}};


    const Shape shape1(points1, "xs1");
    const Shape shape2(points2, "xs2");
    const Shape shape3(points3, "xs3");
    const std::vector shapes = {shape1, shape2, shape3};
    Geometry geometry(shapes, xs_map);

    CHECK(geometry.get_index(-2.0,0,0)==back_index); //outside of all shapes
    CHECK(geometry.get_index(0,0,0)==index3); //should be in the third shape
    CHECK(geometry.get_index(-0.5,0.4,0)==index1);
    CHECK(geometry.get_index(0.6,0.2,0)==index2);

    //overlapping here, but shape3 takes precedence because it was added later
    CHECK(geometry.get_index(0,0.16,0)==index3);
    CHECK(geometry.get_index(0.37,0.2,0)==index3);
    CHECK(geometry.get_index(0.37,0.2,1.4)==back_index);
    CHECK(geometry.get_index(0.37,0.2,1.05)==back_index);
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


TEST_CASE("Check Gy operator") {
    int numy = 1000;
    int numz = 120;
    std::complex<double> amplitude = {0.3, 0.5}; //u0
    std::complex<double> preFactor = {-0.1, 0.8}; //p
    auto ygrid = AuxiliaryFunctions::linspace(-4.0, 4.0, numy);
    auto zgrid = AuxiliaryFunctions::linspace(-2.0, 6.0, numz);
    double k0 = 1.0;
    double reference_index = 1.6;
    double index = 1.5; //using a uniform medium.
    PML pmly(1.0, 5.0, ygrid.front(), ygrid.back());

    //Using some example field: u = u0*sin(k0*y)*sin(2*k0*z)
    multi_array<std::complex<double>, 2> field(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_derivative(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            field[i][j] = amplitude * sin(k0 * ygrid[i]) * sin(2.0 * k0 * zgrid[j]);
            field_derivative[i][j] = amplitude * k0 * cos(k0 * ygrid[i]) * sin(2.0 * k0 * zgrid[j]);
        }
    }

    // (1+p*Gy)u = (1+p/2*k0^2*(n^2-n0^2))*u+p*eta*eta'*u'+p*eta^2*u''
    // where prime is the derivative w.r.t. y, and eta is the PML factor.
    multi_array<std::complex<double>, 2> reference_values(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                reference_values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                auto eta = pmly.get_pml_factor(ygrid[i], index);
                auto eta_prime = pmly.get_pml_factor_derivative(ygrid[i], index);
                auto coef1 = 1.0 + 0.5 * preFactor * k0 * k0 * (index * index - reference_index * reference_index) -
                             preFactor * k0 * k0 * eta * eta;
                auto coef2 = preFactor * eta * eta_prime;
                reference_values[i][j] = coef1 * field[i][j] + coef2 * field_derivative[i][j];
            }
        }
    }

    //now apply the discretized form of the operator.
    multi_array<std::complex<double>, 2> values(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // z value on the mid-point and neighbors is the same, as we only apply the derivative in y direction here.
                values[i][j] = OperatorSuite::apply_right_hand_operator(ygrid[i], ygrid[i - 1], ygrid[i + 1], index,
                                                                        index, index, reference_index, k0, preFactor,
                                                                        field[i][j], field[i - 1][j], field[i + 1][j],
                                                                        &pmly);
            }
        }
    }

    double rms_error = 0.0;
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            double base = std::abs(values[i][j] - reference_values[i][j]);
            rms_error += base * base;
        }
    }
    rms_error = std::sqrt(rms_error / (numy * numz));

    //it will never be perfectly equal anyway, because it's an approximation on a grid. But the error should be small.
    double abs_tol = 1.0e-3;
    CHECK(rms_error <abs_tol);
}

TEST_CASE("Check Gz operator") {
    int numy = 120;
    int numz = 1000;
    std::complex<double> amplitude = {0.3, 0.5}; //u0
    std::complex<double> preFactor = {-0.1, 0.8}; //p
    auto ygrid = AuxiliaryFunctions::linspace(-4.0, 4.0, numy);
    auto zgrid = AuxiliaryFunctions::linspace(-2.0, 6.0, numz);
    double k0 = 1.0;
    double reference_index = 1.6;
    double index = 1.5; //using a uniform medium.
    PML pmlz(1.0, 5.0, zgrid.front(), zgrid.back());

    //Using some example field: u = u0*sin(k0*y)*sin(2*k0*z)
    multi_array<std::complex<double>, 2> field(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_derivative(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            field[i][j] = amplitude * sin(k0 * ygrid[i]) * sin(2.0 * k0 * zgrid[j]);
            field_derivative[i][j] = 2.0 * amplitude * k0 * sin(k0 * ygrid[i]) * cos(2.0 * k0 * zgrid[j]);
        }
    }

    // (1+p*Gz)u = (1+p/2*k0^2*(n^2-n0^2))*u+p*eta*eta'*u'+p*eta^2*u''
    // where prime is the derivative w.r.t. z, and eta is the PML factor.
    // and u'' = -4*k0^2*u
    multi_array<std::complex<double>, 2> reference_values(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                reference_values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                auto eta = pmlz.get_pml_factor(zgrid[j], index);
                auto eta_prime = pmlz.get_pml_factor_derivative(zgrid[j], index);
                auto coef1 = 1.0 + 0.5 * preFactor * k0 * k0 * (index * index - reference_index * reference_index) -
                             4.0 * preFactor * k0 * k0 * eta * eta;
                auto coef2 = preFactor * eta * eta_prime;
                reference_values[i][j] = coef1 * field[i][j] + coef2 * field_derivative[i][j];
            }
        }
    }

    //now apply the discretized form of the operator.
    multi_array<std::complex<double>, 2> values(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // y value on the mid-point and neighbors is the same, as we only apply the derivative in z direction here.
                values[i][j] = OperatorSuite::apply_right_hand_operator(zgrid[j], zgrid[j - 1], zgrid[j + 1], index,
                                                                        index, index, reference_index, k0, preFactor,
                                                                        field[i][j], field[i][j - 1], field[i][j + 1],
                                                                        &pmlz);
            }
        }
    }


    double rms_error = 0.0;
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            double base = std::abs(values[i][j] - reference_values[i][j]);
            rms_error += base * base;
        }
    }
    rms_error = std::sqrt(rms_error / (numy * numz));

    double abs_tol = 1.0e-3;
    CHECK(rms_error <abs_tol);
}

//(1+pGy)(1+pGz)u
TEST_CASE("Check RHS operator") {
    int numy = 1200;
    int numz = 1000;
    std::complex<double> amplitude = {0.3, 0.5}; //u0
    std::complex<double> preFactor = {-0.1, 0.8}; //p
    auto ygrid = AuxiliaryFunctions::linspace(-4.0, 4.0, numy);
    auto zgrid = AuxiliaryFunctions::linspace(-2.0, 6.0, numz);
    double k0 = 1.0;
    double reference_index = 1.6;
    double index = 1.5; //using a uniform medium.
    PML pmly(1.0, 5.0, ygrid.front(), ygrid.back());
    PML pmlz(1.0, 5.0, zgrid.front(), zgrid.back());

    //Using some example field: u = u0*sin(k0*y)*sin(2*k0*z)
    multi_array<std::complex<double>, 2> field(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_derivative_y(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_derivative_z(extents[numy][numz]);
    multi_array<std::complex<double>, 2> field_derivative_yz(extents[numy][numz]);
    multi_array<double, 2> index_array(extents[numy][numz]);
    std::fill_n(index_array.data(), index_array.num_elements(), index);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            field[i][j] = amplitude * sin(k0 * ygrid[i]) * sin(2.0 * k0 * zgrid[j]);
            field_derivative_y[i][j] = amplitude * k0 * cos(k0 * ygrid[i]) * sin(2.0 * k0 * zgrid[j]);
            field_derivative_z[i][j] = 2.0 * amplitude * k0 * sin(k0 * ygrid[i]) * cos(2.0 * k0 * zgrid[j]);
            field_derivative_yz[i][j] = 2.0 * amplitude * k0 * k0 * cos(k0 * ygrid[i]) * cos(2.0 * k0 * zgrid[j]);
        }
    }

    /* RHS = (1+pGy)(1+pGz)u =
     * t^2 u + t*p*xi^2 d^u/dz^2 + t*p*xi*xi'*du/dz +
     * t*p*eta^2 d^2u/dy^2 + p^2*eta^2*xi^2* d^4 u/(dy^2*dz^2) + p^2 * eta^2*xi* xi'*d^3 u/(dy^2*dz) +
     * t*p*eta*eta'*du/dy + p^2*eta*eta'*xi^2* d^3 u/(dy*dz^2) + p^2 * eta*eta'*xi*xi'* d^2 u/(dy*dz).
     * with eta the PML factor in y direction, xi the PML factor in z direction.
     * And t = 1+p/2*k0^2*(n^2-n0^2).
     */
    auto tFactor = 1.0 + 0.5 * preFactor * k0 * k0 * (index * index - reference_index * reference_index);
    multi_array<std::complex<double>, 2> reference_values(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                reference_values[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                auto eta = pmly.get_pml_factor(ygrid[i], index);
                auto eta_prime = pmly.get_pml_factor_derivative(ygrid[i], index);
                auto xi = pmlz.get_pml_factor(zgrid[j], index);
                auto xi_prime = pmlz.get_pml_factor_derivative(zgrid[j], index);
                auto coef1 = tFactor * tFactor - tFactor * preFactor * k0 * k0 * (eta * eta + 4.0 * xi * xi)
                             + 4.0 * preFactor * preFactor * eta * eta * xi * xi * k0 * k0 * k0 * k0;
                auto coef2 = preFactor * eta * eta_prime * (tFactor - 4.0 * k0 * k0 * preFactor * xi * xi);
                auto coef3 = preFactor * xi * xi_prime * (tFactor - k0 * k0 * preFactor * eta * eta);
                auto coef4 = preFactor * preFactor * eta * eta_prime * xi * xi_prime;
                reference_values[i][j] = coef1 * field[i][j] + coef2 * field_derivative_y[i][j] + coef3 *
                                         field_derivative_z[i][j] + coef4 * field_derivative_yz[i][j];
            }
        }
    }

    //now apply the discretized form of the operator.
    auto values = OperatorSuite::get_rhs(field, ygrid, zgrid, index_array, reference_index, k0, preFactor, &pmly,
                                         &pmlz);

    double rms_error = 0.0;
    // skip first and last two rows/columns of points. The first because it is zero anyway,
    // the second because you get a very large numerical gradient there. The solution snaps to zero there.
    // the reference case does not account for that because the derivative is taken analytically, so only justified to compare in the bulk.
    // Note, this spike only appears on the edge of the domain in the y direction. This is because the rhs is calculated in two stages,
    // In the first step the edges of the field were not zerod, so the numerical derivative is smooth. Then the edges get zeroed
    // and the derivative in y direction is taken, resulting in the spike. This does not appear in the previous tests of Gy, Gz, because there the
    // field is not zeroed on the edge. Note that this spike is not a numerical issue, it would also happen on the left hand side of the equation so it is in balance.
    for (int i = 2; i < numy - 2; i++) {
        for (int j = 2; j < numz - 2; j++) {
            double base = std::abs(values[i][j] - reference_values[i][j]);
            rms_error += base * base;
        }
    }
    rms_error = std::sqrt(rms_error / (numy * numz));

    double abs_tol = 1.0e-3;
    CHECK(rms_error <abs_tol);
}
