#include <cmath>
#include <algorithm>
#include "PCSEngine.hpp"
#include "FftComplex.hpp"
#include "G2INIT.h"

PCSEngine::PCSEngine(double sampleFreq)
    : m_sampleFreq(sampleFreq), 
      m_nco(10, (float)sampleFreq), // Build trig tables once
      m_workspace(16384),           // Pre-allocate buffers
      m_accumulatedMag(16384) 
{}

void PCSEngine::initPrn(int prn)
{
    if (codeFfts.find(prn) != codeFfts.end())
        return;

    std::vector<std::complex<double>> codeVec(N);
    G2INIT sv(prn, 0);
    for (size_t idx = 0; idx < N; idx++)
    {
        int chipIdx = (idx / 16) % 1023;
        double codeVal = (double)sv.CODE[chipIdx];
        codeVec[idx] = std::complex<double>(codeVal, 0.0);
    }
    Fft::transform(codeVec, false);
    for (auto &val : codeVec)
        val = std::conj(val);
    
    codeFfts[prn] = std::move(codeVec);
}

AcqResult PCSEngine::search(int prn,
                            const std::vector<std::complex<double>> &rawData,
                            double centerFreq, int binRange, float binWidth)
{
    AcqResult bestResult = {0, 0, -99.0, -99.0};
    int numBlocks = (int)(rawData.size() / N);
    if (numBlocks == 0)
        return bestResult;

    initPrn(prn);
    
    // Cache the PRN FFT reference to avoid map lookups in the loop
    const std::vector<std::complex<double>>& currentCodeFft = codeFfts[prn];

    for (int bin = -binRange; bin <= binRange; bin++)
    {
        // Reuse pre-allocated buffer instead of declaring std::vector here
        std::fill(m_accumulatedMag.begin(), m_accumulatedMag.end(), 0.0);
        
        m_nco.SetFrequency((float)(centerFreq + (bin * binWidth)));
        m_nco.m_phase = 0; // Reset phase for consistent search start

        for (int b = 0; b < numBlocks; b++)
        {
            const std::complex<double>* rawPtr = &rawData[b * N];

            for (size_t idx = 0; idx < N; idx++)
            {
                uint32_t ncoIdx = m_nco.clk();
                // Avoid creating temporary std::complex objects where possible
                double c = (double)m_nco.cosine(ncoIdx);
                double s = (double)m_nco.sine(ncoIdx);
                
                // Manual complex multiply: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
                double re = rawPtr[idx].real();
                double im = rawPtr[idx].imag();
                
                m_workspace[idx] = { re * c - im * s, re * s + im * c };
            }

            Fft::transform(m_workspace, false);
            
            for (size_t idx = 0; idx < N; idx++)
                m_workspace[idx] *= currentCodeFft[idx];
            
            Fft::transform(m_workspace, true);

            for (size_t idx = 0; idx < N; idx++)
            {
                // std::abs() is the magnitude; divide by N for FFT normalization
                m_accumulatedMag[idx] += std::abs(m_workspace[idx]) / N;
            }
        }

        double maxMag = 0;
        int peakIndex = 0;
        double sumPower = 0;
        for (size_t idx = 0; idx < N; idx++)
        {
            double mag = m_accumulatedMag[idx] / numBlocks;
            double pwr = mag * mag;
            sumPower += pwr;
            if (mag > maxMag)
            {
                maxMag = mag;
                peakIndex = (int)idx;
            }
        }

        double peakPower = maxMag * maxMag;
        double avgNoisePower = (sumPower - peakPower) / (double)(N - 1);
        
        // Protect against log of zero/negative
        double snr = (avgNoisePower > 0) ? 10.0 * std::log10(peakPower / avgNoisePower) : -99.0;

        if (snr > bestResult.snr)
        {
            bestResult = {bin, peakIndex, maxMag, snr};
        }
    }
    return bestResult;
}