#include <iostream>
#include "Engine.h"

int main(const int argc, char *argv[]) {
    if (argc == 1) {
        std::cout << "Insufficient command line arguments, json input file name missing!" << std::endl;
        return 1;
    } else if (argc > 2) {
        std::cout << "Too many command line arguments!" << std::endl;
        return 2;
    } else {
        const std::string inputFileName(argv[1]);
        Engine engine(inputFileName); //initialize engine
        engine.run(); //start the solving process
        return 0;
    }
}