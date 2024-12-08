//
// Created by mike on 11/17/24.
//

#include "ProgressBar.h"
#include <iostream>
#include <format>

ProgressBar::ProgressBar(const int max_steps) : max_steps(max_steps), checkpoint_1percent(max_steps / 100),
                                                checkpoint_10percent(max_steps / 10),
                                                checkpoint_50percent(max_steps / 2) {
    begin = std::chrono::steady_clock::now();
}

void ProgressBar::update(const int current_step) const {
    if (current_step == checkpoint_1percent) {
        print_progress(1.0);
    }
    if (current_step == checkpoint_10percent) {
        print_progress(10.0);
    }
    if (current_step == checkpoint_50percent) {
        print_progress(50.0);
    }
}

void ProgressBar::finalize() const {
    print_progress(100.0);
}

void ProgressBar::print_progress(double percentage) const {
    const auto end = std::chrono::steady_clock::now();
    auto delta = 1.0e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).
                     count());
    std::cout << std::format("{0:d}% reached,\t elapsed time = {1:.3f} (s),\t expected total run time = {2:.3f} (s).\n",
                             static_cast<int>(percentage), delta, delta * 100.0 / percentage);
}
