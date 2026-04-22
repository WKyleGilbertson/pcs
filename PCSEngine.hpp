#ifndef PCSEngine_HPP
#define PCSEngine_HPP

#include <vector>
#include <complex>
#include <map>
#include "NCO.h" // Must include this for m_nco member

struct AcqResult {
    int bin;
    int peakIndex;
    double peakMagnitude;
    double snr;
};

class PCSEngine {
private:
    const size_t N = 16384;
    double m_sampleFreq;

    // These MUST be here for the .cpp to find them:
    NCO m_nco; 
    std::vector<std::complex<double>> m_workspace;
    std::vector<double> m_accumulatedMag;

    std::map<int, std::vector<std::complex<double>>> codeFfts;

public:
    PCSEngine(double sampleFreq);
    void initPrn(int prn);
    AcqResult search(int prn, const std::vector<std::complex<double>> &rawData, 
                     double centerFreq, int binRange, float binWidth);
};

#endif