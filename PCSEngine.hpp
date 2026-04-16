#ifndef PCSEngine_HPP
#define PCSEngine_HPP

#include <vector>
#include <complex>
#include <map>

struct AcqResult
{
    int bin;
    int peakIndex;
    double peakMagnitude;
    double snr;
};

class PCSEngine
{
private:
    const size_t N = 16384;
    double m_sampleFreq;
    std::map<int, std::vector<std::complex<double>>> codeFfts;

public:
    PCSEngine(double sampleFreq);

    void initPrn(int prn);

    AcqResult search(int prn,
                     const std::vector<std::complex<double>> &rawData,
                     double centerFreq,
                     int binRange,
                     float binWidth);
};

#endif
