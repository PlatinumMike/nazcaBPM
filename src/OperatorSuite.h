//
// Created by mike on 10/15/24.
//

#ifndef OPERATORSUITE_H
#define OPERATORSUITE_H
#include <complex>
#include "PML.h"
#include <boost/multi_array.hpp>
using boost::multi_array;


class OperatorSuite {
public:
    /**
     * applies the operator (1+prefactorGy) or (1+prefactorGz). You can decide which one to do based on the neighbor values you pass into the function.
     * @param position_mid mid-point
     * @param position_back position of one grid point back in the direction of the derivative.
     * @param position_forward position of one grid point forward in the direction of the derivative.
     * @param index_mid refractive index at the mid-point.
     * @param index_back refractive index one grid point back.
     * @param index_forward refractive index one grid point forward.
     * @param reference_index reference refractive index.
     * @param k0 vacuum wavenumber.
     * @param preFactor complex coefficient to multiply the operator with.
     * @param fieldValue_mid field value at the mid-point
     * @param fieldValue_back field value of one grid point back in the direction of the derivative.
     * @param fieldValue_forward field value of one grid point forward in the direction of the derivative.
     * @param pmlPtr point to PML object
     * @return operator value at the given point of interest
     */
    static std::complex<double> apply_right_hand_operator(double position_mid, double position_back,
                                                          double position_forward,
                                                          double index_mid,
                                                          double index_back,
                                                          double index_forward, double reference_index,
                                                          double k0,
                                                          std::complex<double> preFactor,
                                                          std::complex<double> fieldValue_mid,
                                                          std::complex<double> fieldValue_back,
                                                          std::complex<double> fieldValue_forward, const PML *pmlPtr);

    static multi_array<std::complex<double>, 2> get_rhs(const multi_array<std::complex<double>, 2> &field,
                                                        const std::vector<double> &ygrid,
                                                        const std::vector<double> &zgrid,
                                                        const multi_array<double, 2> &index,
                                                        double reference_index, double k0,
                                                        std::complex<double> preFactor,
                                                        const PML *pmlyPtr, const PML *pmlzPtr);


    static std::vector<std::complex<double> > solve_system(const std::vector<double> &position_mid,
                                                           const std::vector<double> &position_back,
                                                           const std::vector<double> &position_forward,
                                                           const std::vector<double> &index_mid,
                                                           const std::vector<double> &index_back,
                                                           const std::vector<double> &index_forward,
                                                           double reference_index,
                                                           double k0,
                                                           std::complex<double> preFactor,
                                                           const std::vector<std::complex<double> > &rhs_slice,
                                                           const PML *pmlPtr);
};


#endif //OPERATORSUITE_H
