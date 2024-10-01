//
// Created by mike on 10/1/24.
//

#ifndef READERS_H
#define READERS_H

#include "../Parameters.h"
#include <string>


class Readers{
public:
    static Parameters readJSON(const std::string &inputFileName);

private:
    /**
     * Do some basic tests on the inputs before launching a simulation.
     * This avoids wasting large amounts of compute time on invalid inputs.
     * @param params input parameters
     */
    static void testParameters(const Parameters & params);

};



#endif //READERS_H
