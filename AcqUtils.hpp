#ifndef ACQ_UTILS_HPP
#define ACQ_UTILS_HPP

#include <vector>
#include <string>
//#include "version.h"    // This provides the SWV struct definition
#include "kiss_fft.h"
#include "PCSEngine.hpp"

namespace AcqUtils {

enum class DataFormat { BIN, RAW_USB, RAW_BVF };

    struct Config {
        std::string filename = "IF.bin";
        int numMs = 1;
        std::vector<int> prnsToSearch;
        DataFormat format = DataFormat::BIN;
        bool showHelp = false; // <-- NEW
    };

    void PrintUsage(const char* progName);

    void PrintVersion();
    Config ParseArgs(int argc, char *argv[]);
    bool LoadRawData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs, bool isFNLN);
    bool LoadBinData(const std::string& filename, std::vector<kiss_fft_cpx>& data, int numMs);
    void PrintHeader();
    void PrintResult(int prn, const AcqResult& result);
}
#endif