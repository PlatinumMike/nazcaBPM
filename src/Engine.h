//
// Created by mike on 9/29/24.
//

#ifndef ENGINE_H
#define ENGINE_H

#include <string>

#include "Parameters.h"

/**
 * Note, the coordinate system is described in the latex docs.
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
