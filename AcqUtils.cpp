#include "AcqUtils.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>
namespace AcqUtils {

void PrintUsage(const char* progName) {
    printf("Parallel Code Search (PCS) Tool\n");
    printf("Usage: %s [filename] [options]\n\n", progName);
    printf("Options:\n");
    printf("  -h, -?       Show this help message and exit.\n");
    printf("  -m <ms>      Set coherent integration time in milliseconds (default: 1).\n");
    printf("  -p <prns>    Comma-separated list of PRNs to search (default: 1-32, 131, 133, 135).\n");
    printf("  -f <file>    Explicitly specify input file (positional argument preferred).\n\n");
    printf("Format Options (Default is BIN):\n");
    printf("  -u           Read input as RAW USB data.\n");
    printf("  -b, --bvf    Read input as RAW BeagleV-Fire data.\n\n");
    printf("Example:\n");
    printf("  %s BVF-cap.raw -m 5 -b\n", progName);
}

bool LoadRawData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs, bool isBVF) {
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
        
auto map_bits = [](uint8_t mag, uint8_t sign) -> int16_t {
            return (mag ? 3 : 1) * (sign ? -1 : 1);
        };

        for (size_t i = 0; i < bytesToRead; i++) {
            uint8_t b = ingest_buf[i];
            kiss_fft_cpx& s0 = data[offset + (2 * i)];
            kiss_fft_cpx& s1 = data[offset + (2 * i) + 1];

            if (isBVF) {
                // BVF Layout: {i1_s, i1_m, q1_s, q1_m, i0_s, i0_m, q0_s, q0_m}
                s0.i = map_bits((b >> 0) & 1, (b >> 1) & 1); // Q (Imaginary)
                s0.r = map_bits((b >> 2) & 1, (b >> 3) & 1); // I (Real)
                s1.i = map_bits((b >> 4) & 1, (b >> 5) & 1); // Q (Imaginary)
                s1.r = map_bits((b >> 6) & 1, (b >> 7) & 1); // I (Real)
            } else {
                // USB Layout (isBVF == false)
                unpackL1IF(b, s0, s1, true); 
            }
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

    // Trigger help if no arguments are provided
    if (argc == 1) {
        cfg.showHelp = true;
        return cfg;
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Handle Help Flags
        if (arg == "-h" || arg == "-?" || arg == "--help") {
            cfg.showHelp = true;
            return cfg;
        }
        // Handle File Flag (Legacy -f support)
        else if (arg.substr(0, 2) == "-f") {
            cfg.filename = (arg.length() > 2) ? arg.substr(2) : argv[++i];
        } 
        // Handle Milliseconds Flag
        else if (arg.substr(0, 2) == "-m") {
            std::string val = (arg.length() > 2) ? arg.substr(2) : argv[++i];
            cfg.numMs = std::stoi(val); 
        }
        // Handle Format Flags
        else if (arg == "-u") {
            cfg.format = DataFormat::RAW_USB;
        }
        else if (arg == "-b" || arg == "--bvf") {
            cfg.format = DataFormat::RAW_BVF;
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
        // Handle Positional Filename (Must be last so it catches anything without a dash)
        else if (arg[0] != '-') {
            cfg.filename = arg;
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