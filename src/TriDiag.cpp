//
// Created by mike on 10/1/24.
//

#include "TriDiag.h"

std::vector<double> TriDiag::solve_thomas(const std::vector<double> &lower_diag, std::vector<double> &diag,
                                          const std::vector<double> &upper_diag, std::vector<double> &rhs) {
    const size_t size = diag.size();
    std::vector<double> solution(size);

    //forward sweep
    double w;
    for (size_t i = 0; i < size - 1; i++) {
        w = lower_diag[i] / diag[i];
        diag[i + 1] -= w * upper_diag[i];
        rhs[i + 1] -= w * rhs[i];
    }

    // back substitution
    solution[size - 1] = rhs[size - 1] / diag[size - 1];
    for (size_t i = 1; i < size; i++) {
        solution[size - 1 - i] = (rhs[size - 1 - i] - upper_diag[size - 1 - i] * solution[size - i]) / diag[
                                     size - 1 - i];
    }

    return solution;
}
