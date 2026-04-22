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

AcqResult PCSEngine::search(int prn, const std::vector<std::complex<double>> &rawData,
                            double centerFreq, int binRange, float binWidth)
{
    AcqResult bestResult = {0, 0, -99.0, -99.0};
    int numBlocks = (int)(rawData.size() / N);
    if (numBlocks == 0)
        return bestResult;

    initPrn(prn);
    const std::vector<kiss_fft_cpx> &currentCodeFft = codeFfts[prn];

    for (int bin = -binRange; bin <= binRange; bin++)
    {
        std::fill(m_accumulatedMag.begin(), m_accumulatedMag.end(), 0.0);

        // Update NCO for this Doppler bin
        m_nco.SetFrequency((float)(centerFreq + (bin * binWidth)));
        m_nco.m_phase = 0;

        for (int b = 0; b < numBlocks; b++)
        {
            // Safety check: Pointer to the start of this block
            const std::complex<double> *blockStart = &rawData[b * N];

            for (size_t idx = 0; idx < N; idx++)
            {
                uint32_t ncoIdx = m_nco.clk();
                double c = (double)m_nco.cosine(ncoIdx);
                double s = (double)m_nco.sine(ncoIdx);

                // Explicitly pull real/imag to avoid any hidden padding issues
                double re = blockStart[idx].real();
                double im = blockStart[idx].imag();

                m_workspace[idx].r = re * c - im * s;
                m_workspace[idx].i = re * s + im * c;
            }

            // Forward FFT (Time -> Freq)
            kiss_fft(m_cfg_fwd, m_workspace.data(), m_workspace.data());

            // Point-wise complex multiplication (Correlation in Freq Domain)
            for (size_t idx = 0; idx < N; idx++)
            {
                double a = m_workspace[idx].r;
                double b = m_workspace[idx].i; // Fixed from .imag to .i
                double c = currentCodeFft[idx].r;
                double d = currentCodeFft[idx].i;

                // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
                m_workspace[idx].r = a * c - b * d;
                m_workspace[idx].i = a * d + b * c;
            }

            // Inverse FFT (Freq -> Time correlation)
            kiss_fft(m_cfg_inv, m_workspace.data(), m_workspace.data());

            for (size_t idx = 0; idx < N; idx++)
            {
                // Magnitude squared is faster and avoids sqrt in the inner loop
                double r = m_workspace[idx].r;
                double i = m_workspace[idx].i;
                double mag = std::sqrt(r * r + i * i);
                m_accumulatedMag[idx] += mag;
            }
        }

        // Peak detection (mag is currently sum of mags)
        double maxMag = 0;
        int peakIndex = 0;
        double sumPower = 0;

        for (size_t idx = 0; idx < N; idx++)
        {
            // Correct normalization for KissFFT:
            // The IFFT result is N times larger than the correlation result.
            double mag = m_accumulatedMag[idx] / (double)(N * numBlocks);
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
        double snr = (avgNoisePower > 1e-18) ? 10.0 * std::log10(peakPower / avgNoisePower) : -99.0;

        if (snr > bestResult.snr)
        {
            bestResult = {bin, peakIndex, maxMag, snr};
        }
    }
    return bestResult;
}