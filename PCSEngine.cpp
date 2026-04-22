#include <cmath>
#include <algorithm>
#include <iostream>
#include "PCSEngine.hpp"
#include "G2INIT.h"

#undef KISS_FFT_ERROR
extern "C"
{
    void KISS_FFT_ERROR(const char *msg)
    {
        std::cerr << "KissFFT Error: " << msg << std::endl;
    }
}

PCSEngine::PCSEngine(double sampleFreq)
    : m_sampleFreq(sampleFreq),
      m_nco(10, (float)sampleFreq),
      m_workspace(16384),
      m_accumulatedMag(16384)
{
    // Allocate KissFFT plans
    m_cfg_fwd = kiss_fft_alloc(16384, 0, NULL, NULL);
    m_cfg_inv = kiss_fft_alloc(16384, 1, NULL, NULL);
}

PCSEngine::~PCSEngine()
{
    free(m_cfg_fwd);
    free(m_cfg_inv);
}

void PCSEngine::initPrn(int prn)
{
    if (codeFfts.find(prn) != codeFfts.end())
        return;

    std::vector<kiss_fft_cpx> codeVec(N);
    G2INIT sv(prn, 0);
    for (size_t idx = 0; idx < N; idx++)
    {
        int chipIdx = (idx / 16) % 1023;
        codeVec[idx].r = (double)sv.CODE[chipIdx];
        codeVec[idx].i = 0.0;
    }

    // Forward FFT for the code
    kiss_fft(m_cfg_fwd, codeVec.data(), codeVec.data());

    // Conjugate for the circular correlation
    for (size_t idx = 0; idx < N; idx++)
    {
        codeVec[idx].i = -codeVec[idx].i;
    }

    codeFfts[prn] = std::move(codeVec);
}

AcqResult PCSEngine::search(int prn, const std::vector<std::complex<float>> &rawData,
                            float centerFreq, int binRange, float binWidth)
{
    AcqResult bestResult = {0, 0, -99.0f, -99.0f};
    int numBlocks = (int)(rawData.size() / N);
    if (numBlocks == 0)
        return bestResult;

    initPrn(prn);
    const std::vector<kiss_fft_cpx> &currentCodeFft = codeFfts[prn];

    for (int bin = -binRange; bin <= binRange; bin++)
    {
        std::fill(m_accumulatedMag.begin(), m_accumulatedMag.end(), 0.0f);

        // Update NCO for this Doppler bin
        m_nco.SetFrequency((float)(centerFreq + (bin * binWidth)));
        m_nco.m_phase = 0;

        for (int b = 0; b < numBlocks; b++)
        {
            // Safety check: Pointer to the start of this block
            const std::complex<float> *blockStart = &rawData[b * N];

            for (size_t idx = 0; idx < N; idx++)
            {
                uint32_t ncoIdx = m_nco.clk();
                float c = (float)m_nco.cosine(ncoIdx);
                float s = (float)m_nco.sine(ncoIdx);

                // Explicitly pull real/imag to avoid any hidden padding issues
                float re = (float)blockStart[idx].real();
                float im = (float)blockStart[idx].imag();

                m_workspace[idx].r = re * c - im * s;
                m_workspace[idx].i = re * s + im * c;
            }

            // Forward FFT (Time -> Freq)
            kiss_fft(m_cfg_fwd, m_workspace.data(), m_workspace.data());

            // Point-wise complex multiplication (Correlation in Freq Domain)
            for (size_t idx = 0; idx < N; idx++)
            {
                float a = m_workspace[idx].r;
                float b = m_workspace[idx].i;
                float c = currentCodeFft[idx].r;
                float d = currentCodeFft[idx].i;

                // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
                m_workspace[idx].r = a * c - b * d;
                m_workspace[idx].i = a * d + b * c;
            }

            // Inverse FFT (Freq -> Time correlation)
            kiss_fft(m_cfg_inv, m_workspace.data(), m_workspace.data());

            for (size_t idx = 0; idx < N; idx++)
            {
                // Magnitude squared is faster and avoids sqrt in the inner loop
                float r = m_workspace[idx].r;
                float i = m_workspace[idx].i;
                float mag = std::sqrt(r * r + i * i);
                m_accumulatedMag[idx] += mag;
            }
        }

        // Peak detection (mag is currently sum of mags)
        float maxMag = 0;
        int peakIndex = 0;
        float sumPower = 0;

        for (size_t idx = 0; idx < N; idx++)
        {
            // Correct normalization for KissFFT:
            // The IFFT result is N times larger than the correlation result.
            float mag = m_accumulatedMag[idx] / (float)(N * numBlocks);
            float pwr = mag * mag;
            sumPower += pwr;

            if (mag > maxMag)
            {
                maxMag = mag;
                peakIndex = (int)idx;
            }
        }

        float peakPower = maxMag * maxMag;
        float avgNoisePower = (sumPower - peakPower) / (float)(N - 1);
        float snr = (avgNoisePower > 1e-18) ? 10.0f * std::log10(peakPower / avgNoisePower) : -99.0f;

        if (snr > bestResult.snr)
        {
            bestResult = {bin, peakIndex, maxMag, snr};
        }
    }
    return bestResult;
}