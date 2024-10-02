//
// Created by mike on 9/30/24.
//

#include "Writers.h"
#include <iostream>
#include <fstream>

void Writers::write_to_csv(std::vector<std::complex<double> > &data, const std::string &filename, int Nx, int Ny) {
    std::cout << "Writing to file " << filename << std::endl;
    //TODO: instead of printing this it would be nice to store this in the file, perhaps with some header info.
    std::cout << "Nx = " << Nx << ", Ny = " << Ny << std::endl;
    std::ofstream file;
    file.open(filename);
    for (const std::complex<double> value: data) {
        file << std::real(value) << ',' << std::imag(value) << '\n';
    }
    file.close();
    std::cout << "Done writing to file " << filename << std::endl;
}

void Writers::write_to_csv(std::vector<double> &data, const std::string &filename, int Nx, int Ny) {
    std::cout << "Writing to file " << filename << std::endl;
    std::cout << "Nx = " << Nx << ", Ny = " << Ny << std::endl;
    std::ofstream file;
    file.open(filename);
    for (const double value: data) {
        file << value << '\n';
    }
    file.close();
    std::cout << "Done writing to file " << filename << std::endl;
}
