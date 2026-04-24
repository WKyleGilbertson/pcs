#include <cmath>
#include <algorithm>
#include <iostream>
#include "PCSEngine.hpp"
#include "G2INIT.h"

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
      m_accumulatedMag(16384),
      m_codeFftCurrent(16384)
{
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

    G2INIT sv(prn, 0);

    // 1. Fill 1ms exactly (16368 samples)
    for (size_t idx = 0; idx < 16368; idx++)
    {
        m_codeFftCurrent[idx].r = (float)sv.CODE[idx / 16];
        m_codeFftCurrent[idx].i = 0.0f;
    }

    // 2. ZERO-PAD the remaining 16 samples to reach N=16384
    // This prevents the "Franken-chip" interference
    for (size_t idx = 16368; idx < N; idx++)
    {
        m_codeFftCurrent[idx].r = 0.0f;
        m_codeFftCurrent[idx].i = 0.0f;
    }

    // 3. Perform Forward FFT
    kiss_fft(m_cfg_fwd, m_codeFftCurrent.data(), m_codeFftCurrent.data());

    // 4. Conjugate (flipped Imaginary)
    for (size_t idx = 0; idx < N; idx++)
    {
        m_codeFftCurrent[idx].i *= -1.0f;
    }

    // 5. Store the clean, padded, frequency-domain replica
    codeFfts[prn] = m_codeFftCurrent;
}

AcqResult PCSEngine::search(int prn, const std::vector<kiss_fft_cpx> &rawData,
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

        m_nco.SetFrequency(centerFreq + (bin * binWidth));

        for (int b = 0; b < numBlocks; b++)
        {
            // Now rawData is already kiss_fft_cpx, so no complex conversion needed
           // const kiss_fft_cpx *blockStart = &rawData[b * 16368];
            const kiss_fft_cpx *blockStart = &rawData[b * N];

            for (size_t idx = 0; idx < 16368; idx++)
            {
                uint32_t ncoIdx = m_nco.clk();
                float c = (float)m_nco.cosine(ncoIdx);
                float s = (float)m_nco.sine(ncoIdx);

                float re = blockStart[idx].r;
                float im = blockStart[idx].i;

                m_workspace[idx].r = re * c - im * s;
                m_workspace[idx].i = re * s + im * c;
            }
            // Zero-pad the remaining 16 samples to reach N=16384
            for (size_t idx = 16368; idx < N; idx++)
            {
                m_workspace[idx].r = 0.0f;
                m_workspace[idx].i = 0.0f;
            }

            kiss_fft(m_cfg_fwd, m_workspace.data(), m_workspace.data());

            for (size_t idx = 0; idx < N; idx++)
            {
                float a = m_workspace[idx].r;
                float b = m_workspace[idx].i;
                float c = currentCodeFft[idx].r;
                float d = currentCodeFft[idx].i;

                m_workspace[idx].r = a * c - b * d;
                m_workspace[idx].i = a * d + b * c;
            }

            kiss_fft(m_cfg_inv, m_workspace.data(), m_workspace.data());

            for (size_t idx = 0; idx < N; idx++)
            {
                float r = m_workspace[idx].r;
                float i = m_workspace[idx].i;
                m_accumulatedMag[idx] += std::sqrtf(r * r + i * i);
            }
        }

        float maxMag = 0.0f;
        int peakIndex = 0;
        float sumPower = 0.0f;
        const float norm = 1.0f / (float)(N * numBlocks);

        for (size_t idx = 0; idx < N; idx++)
        {
            float mag = m_accumulatedMag[idx] * norm;
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
        float snr = (avgNoisePower > 1e-12f) ? 10.0f * std::log10f(peakPower / avgNoisePower) : -99.0f;

        if (snr > bestResult.snr)
        {
            bestResult = {bin, peakIndex, maxMag, snr};
        }
    }
    return bestResult;
}