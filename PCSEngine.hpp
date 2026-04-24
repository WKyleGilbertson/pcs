#ifndef PCSEngine_HPP
#define PCSEngine_HPP

#include <vector>
#include <complex>
#include <map>
#include <cstdint>
#include "NCO.h" // Must include this for m_nco member

#undef kiss_fft_scalar
#define kiss_fft_scalar float
#include "kiss_fft.h"

struct AcqResult
{
    int bin;
    int peakIndex;
    double peakMagnitude;
    double snr;
};

struct CpxInt16 {
    int16_t r;
    int16_t i;
};

class PCSEngine
{
private:
    double m_sampleFreq;
    NCO m_nco;
    int N = 16384; // 2^14 = 16384, 16 more than 16386... zero padding for FFT

    std::vector<int8_t> m_rawInt8; // Interleaved int8_t raw data buffer    
    std::vector<CpxInt16> CpxncoBuff; // 16-bit Integer NCO

    std::vector<kiss_fft_cpx> m_workspace;
    std::vector<float> m_accumulatedMag;
    std::vector<kiss_fft_cpx> m_codeFftCurrent; // <--- Add this line
    std::vector<kiss_fft_cpx> m_ncoBuffer; // Initialize to size N (16384)

    kiss_fft_cfg m_cfg_fwd;
    kiss_fft_cfg m_cfg_inv;

    std::map<int, std::vector<kiss_fft_cpx>> codeFfts;

    static inline void complex_mix(kiss_fft_cpx* out, const kiss_fft_cpx* a, 
                                   const kiss_fft_cpx* b, size_t count) 
    {
        for (size_t i = 0; i < count; ++i) {
            float r = a[i].r * b[i].r - a[i].i * b[i].i;
            float im = a[i].r * b[i].i + a[i].i * b[i].r;
            out[i].r = r;
            out[i].i = im;
        }
    }

// In PCSEngine.hpp (or .cpp)
static inline void complex_mix_bridge(kiss_fft_cpx* out, 
                                      const kiss_fft_cpx* raw_float, 
                                      const CpxInt16* nco_int, 
                                      size_t count) 
{
    // Verification scale: 1.0 / 32767.0
    // This brings the int16 NCO back to the 1.0 range of your previous float NCO
    const float scale = 1.0f / 32767.0f;

    for (size_t i = 0; i < count; ++i) {
        // Shadow Test Math: Float (Data) * Int (NCO)
        float r  = (raw_float[i].r * nco_int[i].r) - (raw_float[i].i * nco_int[i].i);
        float im = (raw_float[i].r * nco_int[i].i) + (raw_float[i].i * nco_int[i].r);

        out[i].r = r * scale;
        out[i].i = im * scale;
    }
}

public:
    PCSEngine(double sampleFreq);
    ~PCSEngine();
    void initPrn(int prn);
    AcqResult search(int prn, const std::vector<kiss_fft_cpx> &rawData,
                     float centerFreq, int binRange, float binWidth);
};

#endif