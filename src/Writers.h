//
// Created by mike on 9/30/24.
//

#ifndef WRITERS_H
#define WRITERS_H
#include <complex>
#include <vector>


class Writers {
public:
    // write data to csv as 1D array
    static void write_to_csv(std::vector<std::complex<double>> &data, const std::string& filename, int Nx, int Ny);

};



#endif //WRITERS_H
