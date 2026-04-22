#ifndef PCSEngine_HPP
#define PCSEngine_HPP

#include <vector>
#include <complex>
#include <map>
#include "NCO.h" // Must include this for m_nco member

#define kiss_fft_scalar double
#include "kiss_fft.h"

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
    NCO m_nco; 

    kiss_fft_cfg m_cfg_fwd;
    kiss_fft_cfg m_cfg_inv;

    // These MUST be here for the .cpp to find them:
    std::vector<kiss_fft_cpx> m_workspace;
    std::vector<double> m_accumulatedMag;
    std::map<int, std::vector<kiss_fft_cpx>> codeFfts;

public:
    PCSEngine(double sampleFreq);
    ~PCSEngine();
    void initPrn(int prn);
    AcqResult search(int prn, const std::vector<std::complex<double>> &rawData, 
                     double centerFreq, int binRange, float binWidth);
};

#endif