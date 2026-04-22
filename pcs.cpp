#include <iostream>
#include <vector>
#include <chrono>
#include <future>
#include <cstdio>
#include "PCSEngine.hpp"
#include "AcqUtils.hpp"

int main(int argc, char *argv[])
{
    // 1. Setup config and version info
    AcqUtils::Config config = AcqUtils::ParseArgs(argc, argv);
    AcqUtils::PrintVersion();

    // 2. Load IF data
    std::vector<kiss_fft_cpx> data(16384 * config.numMs);
    if (!AcqUtils::LoadRawData(config.filename, data, config.numMs)) {
        fprintf(stderr, "Error opening %s\n", config.filename.c_str());
        return 1;
    }

    printf("Searching %d ms of data from: %s\n", config.numMs, config.filename.c_str());
    AcqUtils::PrintHeader();

    // 3. Parallel Sky Search
    auto skyStart = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<int, std::future<AcqResult>>> futures;

    for (int prn : config.prnsToSearch) {
        auto searchTask = [prn, &data]() -> AcqResult {
            ::PCSEngine threadEngine(16.368e6); 
            return threadEngine.search(prn, data, 4.092e6f, 20, 500.0f);
        };
        futures.push_back({prn, std::async(std::launch::async, searchTask)});
    }

    // 4. Collect results (waits for each thread as needed)
    for (auto &f : futures) {
        AcqUtils::PrintResult(f.first, f.second.get());
    }

    auto skyEnd = std::chrono::high_resolution_clock::now();
    auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(skyEnd - skyStart).count();

    printf("------------------------------------------------------------------\n");
    printf("Total Sky Search Time: %lld ms for %zu satellites\n", 
           totalDuration, config.prnsToSearch.size());

    return 0;
}