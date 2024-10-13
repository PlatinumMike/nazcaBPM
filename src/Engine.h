//
// Created by mike on 9/29/24.
//

#ifndef ENGINE_H
#define ENGINE_H

#include <string>

#include "Parameters.h"

/**
* Note, the coordinate system is the same as in Fig. 2.1 of the book "Beam propagation method for design of optical waveguide devices", G. Pedrola, 2016.
* So z is the propagation direction, x is normal out of the PIC surface, and y lies tangent to the PIC. This forms a right-handed coordinate system.
* It's a bit awkward because TE mode points in y direction, and TM in x. But anyway, let's be consistent with the book.
*/
class Engine {
public:
    /**
     * @param inputFileName Name of the json input filename
     */
    explicit Engine(const std::string &inputFileName);

    /**
     * Run the actual simulation
     */
    void run() const;

private:
    Parameters _inputs;
};


#endif //ENGINE_H
