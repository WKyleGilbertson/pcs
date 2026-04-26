#include "AcqUtils.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>
namespace AcqUtils {

bool LoadRawData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs, bool isFNLN) {
    FILE *IN = fopen(filename.c_str(), "rb");
    if (!IN) return false;

    const size_t samplesPerMs = 16368;
    const size_t fftSize = 16384;
    const size_t bytesToRead = 8184; 

    // Move buffer to the heap to prevent stack overflow/corruption
    std::vector<uint8_t> ingest_buf(bytesToRead);
    
    for (int ms = 0; ms < numMs; ms++) {
        // fread into the heap-allocated vector's data pointer
        size_t bytesRead = fread(ingest_buf.data(), 1, bytesToRead, IN);
        
        if (bytesRead < bytesToRead) break;

        size_t offset = (size_t)ms * fftSize;
        
        for (size_t i = 0; i < bytesToRead; i++) {
            // Unpack directly from the heap buffer
            unpackL1IF(ingest_buf[i], data[offset + (2 * i)], data[offset + (2 * i) + 1], isFNLN);
        }
        
        // Zero-pad the remaining 16 samples
        for (size_t i = samplesPerMs; i < fftSize; i++) {
            data[offset + i].r = 0;
            data[offset + i].i = 0;
        }
    }
    
    fclose(IN);
    return true;
}

bool LoadBinData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs) {
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
            data[offset + i].r = (int16_t)raw_buffer[2 * i] << 3;
            data[offset + i].i = (int16_t)raw_buffer[2 * i + 1] << 3;
        }
        
        // Pad the remaining samples to reach the 16384 FFT size
        for (size_t i = 16368; i < 16384; i++) {
            data[offset + i].r = 0;
            data[offset + i].i = 0;
        }
    }
    fclose(IN);
    return true;
}

#include <sstream>

Config ParseArgs(int argc, char *argv[]) {
    Config cfg;
    bool prnsSpecified = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Handle File Flag
        if (arg.substr(0, 2) == "-f") {
            cfg.filename = (arg.length() > 2) ? arg.substr(2) : argv[++i];
        } 
        // Handle Milliseconds Flag (Fixed variable name and logic)
        else if (arg.substr(0, 2) == "-m") {
            std::string val = (arg.length() > 2) ? arg.substr(2) : argv[++i];
            cfg.numMs = std::stoi(val); 
        }
        // Handle PRN Flag (Handles -p131 and -p 131,135)
        else if (arg.substr(0, 2) == "-p") {
            prnsSpecified = true;
            std::string prnList = (arg.length() > 2) ? arg.substr(2) : argv[++i];
            std::stringstream ss(prnList);
            std::string segment;

            while (std::getline(ss, segment, ',')) {
                if (!segment.empty()) {
                    cfg.prnsToSearch.push_back(std::stoi(segment));
                }
            }
        }
    }

    // Default Sky Search if no PRNs specified
    if (!prnsSpecified) {
        for (int k = 1; k <= 32; k++) cfg.prnsToSearch.push_back(k);
        cfg.prnsToSearch.push_back(131); 
        cfg.prnsToSearch.push_back(133); 
        cfg.prnsToSearch.push_back(135);
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