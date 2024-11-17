//
// Created by mike on 11/17/24.
//

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H
#include <chrono>


class ProgressBar {
public:
    explicit ProgressBar(int max_steps);

    /**
     * Print rough indication of simulation progress
     * @param current_step current propagation step
     */
    void update(int current_step) const;

    void finalize() const;

private:
    std::chrono::steady_clock::time_point begin;
    const int max_steps;
    const int checkpoint_1percent;
    const int checkpoint_10percent;
    const int checkpoint_50percent;
};


#endif //PROGRESSBAR_H
