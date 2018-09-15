/*
 * A simple frequency domain high-pass filter class
 */
#ifndef __HPF_H__
#define __HPF_H__

#include <e3e_detection.h>

class HPF
{
  public:
    float cutoff;
    float fs;
    int nfft;

    int f_max;
    std::vector<float> weights;

    HPF(float _cutoff, float _fs, int _nfft) : cutoff(_cutoff), fs(_fs), nfft(_nfft)
    {
      f_max = int(ceilf(2 * cutoff / fs * nfft));

      float step = 1.f / f_max;
      for (int m = 1 ; m < f_max ; m++)
        weights.push_back(sqrtf(m * step));
    }
    void process(e3e_complex *signal)
    {
      for (int m = 1 ; m < f_max ; m++)
        signal[m] *= weights[m];
    }
};

#endif // __HPF_H__
