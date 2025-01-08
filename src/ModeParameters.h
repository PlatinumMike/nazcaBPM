//
// Created by mike on 1/8/25.
//

#ifndef MODEPARAMETERS_H
#define MODEPARAMETERS_H

// specify how much output to write.
enum logging_level {
    ERROR, WARNING, INFO, DEBUG
};

struct ModeParams {
    logging_level level;
    // eps_y, eps_z represent a slight offset of the Gaussian splot so it excites also odd modes, this can be used to find the higher order modes.
    // If you only care about the fundamental mode you can safely set these to 0.
    double eps_y;
    double eps_z;
    // standard deviations of the spot size in micron.
    double std_y;
    double std_z;
    //maximum number of iterations allowed before the mode solver gives up.
    int max_iterations;
    // minimum number of iterations to ensure we get rid of any start-up behavior.
    int min_iterations;
    // absolute tolerance in beta before terminating the mode solve
    double absolute_tolerance;
    //step size in x direction. This is typically the same as for the main BPM simulation, but you can specify a different step size.
    double increment_x;
    // if true, the input value of increment_x is ignored, and taken from the BPM grid.
    bool get_increment_from_bpm;
};

#endif //MODEPARAMETERS_H
