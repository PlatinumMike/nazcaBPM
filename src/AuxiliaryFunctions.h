//
// Created by mike on 9/29/24.
//

#ifndef AUXILIARYFUNCTIONS_H
#define AUXILIARYFUNCTIONS_H

#include <vector>

class AuxiliaryFunctions {
public:
    /**
     * Generate vector with linearly spaced values
     * @param start starting value
     * @param stop final value
     * @param count number of entries to fill
     */
    static std::vector<double> linspace(double start, double stop, int count);

};



#endif //AUXILIARYFUNCTIONS_H
