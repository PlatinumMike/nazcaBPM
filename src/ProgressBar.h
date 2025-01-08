//
// Created by mike on 11/17/24.
//

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H
#include <chrono>
#include <vector>

struct Checkpoint {
    int step;
    double percentage;
};

class ProgressBar {
public:
    ProgressBar(int max_steps, std::vector<double> percentages);

    /**
     * Print rough indication of simulation progress
     * @param current_step current propagation step
     */
    void update(int current_step) const;

    void finalize() const;

private:
    std::chrono::steady_clock::time_point begin;
    const int max_steps;
    std::vector<Checkpoint> checkpoints;

    void print_progress(double percentage) const;
};


#endif //PROGRESSBAR_H
