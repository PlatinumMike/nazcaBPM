#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "../src/TriDiag.h"
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
