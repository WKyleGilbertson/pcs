#include <cmath>
#include <algorithm>
#include "PCSEngine.hpp"
#include "FftComplex.hpp"
#include "G2INIT.h"
#include "NCO.h"

PCSEngine::PCSEngine(double sampleFreq)
    : m_sampleFreq(sampleFreq) {}

void PCSEngine::initPrn(int prn)
{
    // Avoid re-calculating if already in map
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
    codeFfts[prn] = codeVec;
}

AcqResult PCSEngine::search(int prn,
                            const std::vector<std::complex<double>> &rawData,
                            double centerFreq, int binRange, float binWidth)
{
    AcqResult bestResult = {0, 0, -99.0, -99.0};
    NCO nco(10, m_sampleFreq);

    int numBlocks = (int)(rawData.size() / N);
    if (numBlocks == 0)
        return bestResult;

    // Ensure PRN is initialized
    initPrn(prn);

    for (int bin = -binRange; bin <= binRange; bin++)
    {
        std::vector<double> accumulatedMag(N, 0.0);
        nco.SetFrequency(centerFreq + (bin * binWidth));

        for (int b = 0; b < numBlocks; b++)
        {
            std::vector<std::complex<double>> workspace(N);
            for (size_t idx = 0; idx < N; idx++)
            {
                uint32_t ncoIdx = nco.clk();
                std::complex<double> lo(nco.cosine(ncoIdx), nco.sine(ncoIdx));
                workspace[idx] = rawData[b * N + idx] * lo;
            }

            Fft::transform(workspace, false);
            for (size_t idx = 0; idx < N; idx++)
                workspace[idx] *= codeFfts[prn][idx];
            Fft::transform(workspace, true);

            for (size_t idx = 0; idx < N; idx++)
            {
                accumulatedMag[idx] += std::abs(workspace[idx]) / N;
            }
        }

        double maxMag = 0;
        int peakIndex = 0;
        double sumPower = 0;
        for (size_t idx = 0; idx < N; idx++)
        {
            double mag = accumulatedMag[idx] / numBlocks;
            sumPower += (mag * mag);
            if (mag > maxMag)
            {
                maxMag = mag;
                peakIndex = idx;
            }
        }

        double peakPower = maxMag * maxMag;
        double avgNoisePower = (sumPower - peakPower) / (N - 1);
        double snr = 10.0 * log10(peakPower / avgNoisePower);

        if (snr > bestResult.snr)
        {
            bestResult = {bin, peakIndex, maxMag, snr};
        }
    }
    return bestResult;
}
