//
// Created by mike on 10/15/24.
//

#include "OperatorSuite.h"

#include "TriDiag.h"
using boost::extents;

std::complex<double> OperatorSuite::apply_right_hand_operator(const double position_mid, const double position_back,
                                                              const double position_forward,
                                                              const double index_mid,
                                                              const double index_back,
                                                              const double index_forward, const double reference_index,
                                                              const double k0,
                                                              const std::complex<double> preFactor,
                                                              const std::complex<double> fieldValue_mid,
                                                              const std::complex<double> fieldValue_back,
                                                              const std::complex<double> fieldValue_forward,
                                                              const PML *pmlPtr) {
    const double neighbor_distance = position_forward - position_mid;
    auto pmlfactor_mid = pmlPtr->get_pml_factor(position_mid, index_mid);
    auto pmlfactor_back = pmlPtr->get_pml_factor(position_back, index_back);
    auto pmlfactor_forward = pmlPtr->get_pml_factor(position_forward, index_forward);

    std::complex<double> value = fieldValue_mid;
    const std::complex<double> coefficient = preFactor * pmlfactor_mid / (2.0 * neighbor_distance * neighbor_distance);

    // add index contribution
    value += preFactor * 0.5 * k0 * k0 * (index_mid * index_mid - reference_index * reference_index) * fieldValue_mid;
    // add double derivative contribution
    value += coefficient * ((pmlfactor_mid + pmlfactor_forward) * fieldValue_forward - (
                                pmlfactor_forward + 2.0 * pmlfactor_mid + pmlfactor_back) * fieldValue_mid + (
                                pmlfactor_mid + pmlfactor_back) * fieldValue_back);
    return value;
}


multi_array<std::complex<double>, 2> OperatorSuite::get_rhs(const multi_array<std::complex<double>, 2> &field,
                                                            const std::vector<double> &ygrid,
                                                            const std::vector<double> &zgrid,
                                                            const multi_array<double, 2> &index,
                                                            const double reference_index, const double k0,
                                                            const std::complex<double> preFactor,
                                                            const PML *pmlyPtr, const PML *pmlzPtr) {
    int numy = static_cast<int>(ygrid.size());
    int numz = static_cast<int>(zgrid.size());

    //right hand side vector
    multi_array<std::complex<double>, 2> rhs1(extents[numy][numz]);
    // apply the operator (1+p*Gz)
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                rhs1[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // x value on the mid-point and neighbors is the same, as we only apply the derivative in y direction here.
                rhs1[i][j] = apply_right_hand_operator(zgrid[j], zgrid[j - 1],
                                                       zgrid[j + 1], index[i][j], index[i][j - 1], index[i][j + 1],
                                                       reference_index, k0, preFactor,
                                                       field[i][j], field[i][j - 1], field[i][j + 1],
                                                       pmlzPtr);
            }
        }
    }
    //now apply the operator (1+p*Gy)
    multi_array<std::complex<double>, 2> rhs2(extents[numy][numz]);
    for (int i = 0; i < numy; i++) {
        for (int j = 0; j < numz; j++) {
            if (i == 0 || j == 0 || i == numy - 1 || j == numz - 1) {
                //ignore boundary values.
                rhs2[i][j] = std::complex<double>{0.0, 0.0};
            } else {
                // y value on the mid-point and neighbors is the same, as we only apply the derivative in x direction here.
                rhs2[i][j] = apply_right_hand_operator(ygrid[i], ygrid[i - 1], ygrid[i + 1], index[i][j],
                                                       index[i - 1][j], index[i + 1][j], reference_index, k0, preFactor,
                                                       rhs1[i][j], rhs1[i - 1][j], rhs1[i + 1][j], pmlyPtr);
            }
        }
    }
    return rhs2;
}


std::vector<std::complex<double> > OperatorSuite::solve_system(const std::vector<double> &position_mid,
                                                               const std::vector<double> &position_back,
                                                               const std::vector<double> &position_forward,
                                                               const std::vector<double> &index_mid,
                                                               const std::vector<double> &index_back,
                                                               const std::vector<double> &index_forward,
                                                               const double reference_index,
                                                               const double k0,
                                                               const std::complex<double> preFactor,
                                                               const std::vector<std::complex<double> > &rhs_slice,
                                                               const PML *pmlPtr) {
    const auto size = rhs_slice.size();
    const double neighbor_distance = position_forward[0] - position_mid[0];
    // compute matrix entries
    std::vector<std::complex<double> > lower_diag(size - 1);
    std::vector<std::complex<double> > diag(size);
    std::vector<std::complex<double> > upper_diag(size - 1);
    std::vector<std::complex<double> > temp_rhs = rhs_slice;
    for (auto idx = 0; idx < size; idx++) {
        auto pmlfactor_mid = pmlPtr->get_pml_factor(position_mid[idx], index_mid[idx]);
        auto pmlfactor_previous = pmlPtr->get_pml_factor(position_back[idx], index_back[idx]);
        auto pmlfactor_next = pmlPtr->get_pml_factor(position_forward[idx], index_forward[idx]);
        const std::complex<double> coefficient =
                preFactor * pmlfactor_mid / (2.0 * neighbor_distance * neighbor_distance);


        diag[idx] = 1.0 + preFactor * 0.5 * k0 * k0 * (
                        index_mid[idx] * index_mid[idx] - reference_index * reference_index)
                    - coefficient * (pmlfactor_next + 2.0 * pmlfactor_mid + pmlfactor_previous);
        if (idx > 1) {
            lower_diag[idx - 1] = coefficient * (pmlfactor_mid + pmlfactor_previous);
        }
        if (idx < size - 1) {
            upper_diag[idx] = coefficient * (pmlfactor_mid + pmlfactor_next);
        }
    }
    //solve
    auto solution = TriDiag<std::complex<double> >::solve_thomas(lower_diag, diag, upper_diag, temp_rhs);

    return solution;
}
