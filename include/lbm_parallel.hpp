#pragma once

namespace LBMParallel {
void init_threads();
void run_iteration(bool is_last);
void stop_threads();
}; // namespace LBMParallel
