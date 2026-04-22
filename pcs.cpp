#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <future>
#include "version.h"
#include "PCSEngine.hpp"

// Updated to handle multiple blocks and your specific 16368-sample read
void stuffVector(std::vector<kiss_fft_cpx> &vec, FILE *fp, int numMs)
{
  int8_t buffer[32736]; // 16368 samples * 2 bytes
  for (int ms = 0; ms < numMs; ms++)
  {
    size_t bytesRead = fread(buffer, 1, sizeof(buffer), fp);
    if (bytesRead < 32736)
      break;

    size_t offset = ms * 16384;
    for (size_t idx = 0; idx < 16368; idx++)
    {
      vec[offset + idx].r = (float)buffer[2 * idx];
      vec[offset + idx].i = (float)buffer[2 * idx + 1];
    }
    // Zero-pad the remaining 16 samples of this millisecond block
    for (size_t idx = 16368; idx < 16384; idx++)
    {
      vec[offset + idx].r = 0.0f;
      vec[offset + idx].i = 0.0f;
    }
  }
}

int main(int argc, char *argv[])
{
  std::string filename = "IF.bin";
  std::vector<int> prnsToSearch;
  int numMs = 1; // Default to 1ms
  SWV V = {MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION, 0, "N/A", "N/A", "PCS"};

#ifdef CURRENT_HASH
  strncpy(V.GitCI, CURRENT_HASH, 40);
#endif
#ifdef CURRENT_DATE
  strncpy(V.BuildDate, CURRENT_DATE, 16);
#endif
#ifdef CURRENT_NAME
  strncpy(V.Name, CURRENT_NAME, 10);
#endif

  fprintf(stdout, "Parallel Code Search Algorithm\n");
  fprintf(stdout, "%s: GitCI:%s %s v%.1d.%.1d.%.1d\n",
          V.Name, V.GitCI, V.BuildDate,
          V.Major, V.Minor, V.Patch);

  // Manual Argument Parsing
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "-f" && i + 1 < argc)
      filename = argv[++i];
    else if (arg == "-m" && i + 1 < argc)
      numMs = std::stoi(argv[++i]);
    //        else if (arg == "-p" && i + 1 < argc) prnsToSearch.push_back(std::stoi(argv[++i]));
    else if (arg == "-p" && i + 1 < argc)
    {
      std::string list = argv[++i];
      size_t start = 0;
      size_t end = list.find(',');

      while (end != std::string::npos)
      {
        std::string token = list.substr(start, end - start);
        if (!token.empty())
          prnsToSearch.push_back(std::stoi(token));
        start = end + 1;
        end = list.find(',', start);
      }

      // This part captures the last item (or the only item if no comma exists)
      std::string lastToken = list.substr(start);
      if (!lastToken.empty())
      {
        prnsToSearch.push_back(std::stoi(lastToken));
      }
    }
  }

  if (prnsToSearch.empty())
  {
    for (int i = 1; i <= 32; ++i)
      prnsToSearch.push_back(i);
    prnsToSearch.push_back(131);
    prnsToSearch.push_back(133);
    prnsToSearch.push_back(135);
  }

  PCSEngine PCSEngine(16.368e6);
  // We don't use the global PCSEngine for searching anymore to avoid thread conflicts,
  // but it's fine to keep it here for initialization or delete it.
  std::vector<kiss_fft_cpx> data(16384 * numMs);

  FILE *IN = fopen(filename.c_str(), "rb");
  if (!IN)
  {
    fprintf(stderr, "Error opening %s\n", filename.c_str());
    return 1;
  }
  stuffVector(data, IN, numMs);
  fclose(IN);

  printf("Searching %d ms of data from: %s\n", numMs, filename.c_str());
   
  // Start timing the entire sky search
// Start timing the entire sky search
  auto skyStart = std::chrono::high_resolution_clock::now();
// Start timing the entire sky search

  std::vector<std::pair<int, std::future<AcqResult>>> futures;

  for (int prn : prnsToSearch)
  {
    // Define the task separately. 
    // We use [=, &data] to capture prn by value and data by reference.
    auto searchTask = [prn, &data]() -> AcqResult {
        // Use ::PCSEngine to explicitly refer to the CLASS, not a variable
        ::PCSEngine threadEngine(16.368e6);
        return threadEngine.search(prn, data, 4.092e6f, 20, 500.0f);
    };

    // Move the future creation outside the make_pair call to help MSVC deduction
    std::future<AcqResult> f = std::async(std::launch::async, searchTask);
    futures.push_back(std::make_pair(prn, std::move(f)));
  }

  printf("PRN | Bin | Freq Offset | Peak Index | Chip Phase | SNR (dB)\n");
  printf("------------------------------------------------------------------\n");

  for (auto &f : futures)
  {
    int prn = f.first;
    AcqResult result = f.second.get(); 

    printf("%3d | %4d | %10.0f | %10d | %10.2f | %6.2f\n",
           prn, result.bin, (float)result.bin * 500.0f, result.peakIndex,
           (float)result.peakIndex / 16.0f, result.snr);
  }

  auto skyEnd = std::chrono::high_resolution_clock::now();
  auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(skyEnd - skyStart).count();

  printf("------------------------------------------------------------------\n");
  printf("Total Sky Search Time: %lld ms for %zu satellites\n", totalDuration, prnsToSearch.size());
  return 0;
}