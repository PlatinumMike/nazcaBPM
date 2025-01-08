//
// Created by mike on 11/17/24.
//

#include "ProgressBar.h"
#include <iostream>
#include <format>
#include <algorithm>

ProgressBar::ProgressBar(const int max_steps, std::vector<double> percentages) : max_steps(max_steps) {
    begin = std::chrono::steady_clock::now();
    std::sort(percentages.begin(), percentages.end());
    for (const double percentage: percentages) {
        const int step = static_cast<int>(0.01 * percentage * max_steps);
        checkpoints.push_back({step, percentage});
    }
}

void ProgressBar::update(const int current_step) const {
    for (const auto [step, percentage]: checkpoints) {
        if (current_step == step) {
            print_progress(percentage);
        }
    }
}

void ProgressBar::finalize() const {
    print_progress(100.0);
}

void ProgressBar::print_progress(const double percentage) const {
    const auto end = std::chrono::steady_clock::now();
    auto delta = 1.0e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).
                     count());
    std::cout << std::format("{0:d}% reached,\t elapsed time = {1:.3f} (s),\t expected total run time = {2:.3f} (s).",
                             static_cast<int>(percentage), delta, delta * 100.0 / percentage) << std::endl;
}
