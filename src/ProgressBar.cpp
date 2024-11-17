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
        const auto end = std::chrono::steady_clock::now();
        auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
        std::cout << std::format("1% reached, elapsed time = {} (s), expected total run time = {} (s).\n"
                                 , delta, delta * 100);
    }
    if (current_step == checkpoint_10percent) {
        const auto end = std::chrono::steady_clock::now();
        auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
        std::cout << std::format("10% reached, elapsed time = {} (s), expected total run time = {} (s).\n"
                                 , delta, delta * 10);
    }
    if (current_step == checkpoint_50percent) {
        const auto end = std::chrono::steady_clock::now();
        auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
        std::cout << std::format("50% reached, elapsed time = {} (s), expected total run time = {} (s).\n"
                                 , delta, delta * 2);
    }
}

void ProgressBar::finalize() const {
    const auto end = std::chrono::steady_clock::now();
    auto delta = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
    std::cout << std::format("Run completed, elapsed time = {} (s).\n", delta);
}
