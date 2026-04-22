#include "AcqUtils.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>

// This macro is required to handle the '"value"' quoting from your Makefile.
// It ensures the compiler treats the macros as valid string literals.
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

namespace AcqUtils {

void PrintVersion() {
    SWV V;
    V.Major = MAJOR_VERSION;
    V.Minor = MINOR_VERSION;
    V.Patch = PATCH_VERSION;

    // Use STR() to convert Makefile macros to strings safely
    const char* hashStr = STR(CURRENT_HASH);
    const char* dateStr = STR(CURRENT_DATE);
    const char* nameStr = STR(APP_NAME);

    // Parse hex tag from hash string
    sscanf(hashStr, "%x", &V.GitTag);
    
    // Copy into struct buffers defined in version.h
    strncpy(V.GitCI, hashStr, sizeof(V.GitCI) - 1);
    V.GitCI[sizeof(V.GitCI) - 1] = '\0';
    
    strncpy(V.BuildDate, dateStr, sizeof(V.BuildDate) - 1);
    V.BuildDate[sizeof(V.BuildDate) - 1] = '\0';
    
    strncpy(V.Name, nameStr, sizeof(V.Name) - 1);
    V.Name[sizeof(V.Name) - 1] = '\0';

    fprintf(stdout, "Parallel Code Search (PCS) Algorithm\n");
    fprintf(stdout, "%s GitCI:%.7s %s v%d.%d.%d\n",
            V.Name, V.GitCI, V.BuildDate,
            (int)V.Major, (int)V.Minor, (int)V.Patch);
}

bool LoadRawData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs) {
    FILE *IN = fopen(filename.c_str(), "rb");
    if (!IN) return false;

    // Use a unique name to ensure no collision with other headers
    int8_t raw_buffer[32768]; 
    
    for (int ms = 0; ms < numMs; ms++) {
        // We pass the array name raw_buffer which is a pointer to the start
        size_t bytesRead = fread(raw_buffer, 1, 32736, IN);
        if (bytesRead < 32736) break;

        size_t offset = (size_t)ms * 16384;
        for (size_t i = 0; i < 16368; i++) {
            // Indexing is now safe because raw_buffer is explicitly an array
            data[offset + i].r = (float)raw_buffer[2 * i];
            data[offset + i].i = (float)raw_buffer[2 * i + 1];
        }
        
        // Pad the remaining samples to reach the 16384 FFT size
        for (size_t i = 16368; i < 16384; i++) {
            data[offset + i].r = 0.0f;
            data[offset + i].i = 0.0f;
        }
    }
    fclose(IN);
    return true;
}

Config ParseArgs(int argc, char *argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-f" && i + 1 < argc) cfg.filename = argv[++i];
        else if (arg == "-m" && i + 1 < argc) cfg.numMs = std::stoi(argv[++i]);
    }
    if (cfg.prnsToSearch.empty()) {
        for (int i = 1; i <= 32; i++) cfg.prnsToSearch.push_back(i);
        cfg.prnsToSearch.push_back(131); cfg.prnsToSearch.push_back(133); cfg.prnsToSearch.push_back(135);
    }
    return cfg;
}

void PrintHeader() {
    printf("PRN | Bin | Freq Offset | Peak Index | Chip Phase | SNR (dB)\n");
    printf("------------------------------------------------------------------\n");
}

void PrintResult(int prn, const AcqResult& result) {
    printf("%3d | %4d | %10.0f | %10d | %10.2f | %6.2f\n",
           prn, result.bin, (float)result.bin * 500.0f, result.peakIndex,
           (float)result.peakIndex / 16.0f, result.snr);
}

} // namespace AcqUtils