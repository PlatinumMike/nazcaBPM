//
// Created by mike on 10/1/24.
//

#ifndef TRIDIAG_H
#define TRIDIAG_H
#include <vector>


// class for solving tri-diagonal system of equations
class TriDiag {
public:
    /**
     * Solves tridiagonal matrix using the Thomas algorithm.
     * The diagonal, the right-hand size, and solution are size N,
     * the lower and upper diagonal are size N-1.
    */
    // warning, this will modify the input vectors.
    // you can avoid that but that will take more RAM.
    //TODO: use templating to make this work for complex data too
    static std::vector<double> solve_thomas(const std::vector<double> &lower_diag, std::vector<double> &diag,
                                     const std::vector<double> &upper_diag, std::vector<double> &rhs);
};


#endif //TRIDIAG_H
