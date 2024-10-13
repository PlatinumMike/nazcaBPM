//
// Created by mike on 10/1/24.
//

#ifndef TRIDIAG_H
#define TRIDIAG_H
#include <vector>
#include <cassert>

template<typename T>
// class for solving tri-diagonal system of equations
class TriDiag {
public:
    /**
     * Solves tridiagonal matrix using the Thomas algorithm.
     * The diagonal, the right-hand size, and solution are size N,
     * the lower and upper diagonal are size N-1.
    */
    // warning, this will modify the input vectors diag and rhs.
    // you can avoid that but that will take more RAM.
    static std::vector<T> solve_thomas(const std::vector<T> &lower_diag, std::vector<T> &diag,
                                       const std::vector<T> &upper_diag, std::vector<T> &rhs);
};


template<typename T>
std::vector<T> TriDiag<T>::solve_thomas(const std::vector<T> &lower_diag, std::vector<T> &diag,
                                        const std::vector<T> &upper_diag, std::vector<T> &rhs) {
    const size_t size = diag.size();
    assert(("lower diagonal should be one element shorter than the diagonal", lower_diag.size() == size - 1));
    assert(("upper diagonal should be one element shorter than the diagonal", upper_diag.size() == size - 1));
    assert(("right hand side should be the same length as the diagonal", rhs.size() == size));

    std::vector<T> solution(size);

    //forward sweep
    for (size_t i = 0; i < size - 1; i++) {
        T w = lower_diag[i] / diag[i];
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
#endif //TRIDIAG_H
