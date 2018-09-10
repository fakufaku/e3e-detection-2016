#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>

#include "stft.h"
#include "pyramic.h"

#define NUM_SAMPLES 1024
#define NFFT (2 * NUM_SAMPLES)

/*******************/
/* GLOBAL VARS ETC */
/*******************/

STFT engine_in(NUM_SAMPLES, NFFT, 0, 0, PYRAMIC_CHANNELS_IN, STFT_WINDOW_BOTH);
STFT engine_out(NUM_SAMPLES, NFFT, 0, 0, 1, STFT_WINDOW_BOTH);  // single output channel
float buffer_in[NUM_SAMPLES * PYRAMIC_CHANNELS_IN] = {0};
float buffer_out[NUM_SAMPLES] = {0};
float int2float = 1. / (1 << 15);
int16_t float2int = (1 << 15) - 1;

e3e_complex_vector beamformer((NFFT / 2 + 1) * PYRAMIC_CHANNELS_IN);

float delays[] = {
  3.82790651e-05,  4.59348781e-05,  5.35906911e-05,  5.66530164e-05,
  5.81841790e-05,  6.12465042e-05,  6.89023172e-05,  7.65581302e-05,
  7.17732471e-05,  6.02895275e-05,  4.88058080e-05,  4.42123202e-05,
  4.19155763e-05,  3.73220885e-05,  2.58383689e-05,  1.43546494e-05,
  2.39244157e-05,  2.00965092e-05,  1.62686027e-05,  1.47374401e-05,
  1.39718588e-05,  1.24406962e-05,  8.61278965e-06,  4.78488314e-06,
  1.43546494e-05,  2.58383689e-05,  3.73220885e-05,  4.19155763e-05,
  4.42123202e-05,  4.88058080e-05,  6.02895275e-05,  7.17732471e-05,
  2.39244157e-05,  2.00965092e-05,  1.62686027e-05,  1.47374401e-05,
  1.39718588e-05,  1.24406962e-05,  8.61278965e-06,  4.78488314e-06,
  1.01643954e-20,  6.77626358e-21,  6.77626358e-21,  6.77626358e-21,
  6.77626358e-21,  6.77626358e-21,  6.77626358e-21, -0.00000000e+00
};

/*************************/
/* CREATE THE BEAMFORMER */
/*************************/

void compute_beamforming_weights()
{
  const std::complex<float> j(0, 1);  // imaginary number
  float pi_f = (float)M_PI;
  float denom = 1. / PYRAMIC_CHANNELS_IN;

  for (size_t f = 0 ; f < (NFFT / 2 + 1) ; f++)
  {
    float f_hz = (float)f / NFFT * PYRAMIC_SAMPLERATE;
    for (size_t ch = 0 ; ch < PYRAMIC_CHANNELS_IN ; ch++)
      beamformer[f * PYRAMIC_CHANNELS_IN + ch] = std::exp(j * 2.f * pi_f * f_hz * delays[ch]) * denom;
  }
}

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init()
{
  compute_beamforming_weights();

  // Zero the DC and middle element of the output frequency buffer
  engine_out.freq_buffer[0] = 0.;
  engine_out.freq_buffer[NFFT / 2] = 0.;
}

/***********************************/
/* USER-DEFINED PROCESSING ROUTINE */
/***********************************/

// All the processing can be concentrated in this function
void processing(buffer_t &input, buffer_t &output)
{
  for (size_t n = 0 ; n < NUM_SAMPLES * PYRAMIC_CHANNELS_IN ; n++)
    buffer_in[n] = int2float * input[n];

  e3e_complex *spectrum_in = engine_in.analysis(buffer_in);

  // Copy the first two input channels to the output
  for (size_t n = 1 ; n < NFFT / 2 ; n++)
  {
    // The last microphone is not delayed, simply copy
    engine_out.freq_buffer[n] = spectrum_in[PYRAMIC_CHANNELS_IN * n + PYRAMIC_CHANNELS_IN - 1];

    // now apply delay and sum the other channels
    for (size_t ch = 0 ; ch < PYRAMIC_CHANNELS_IN - 1 ; ch++)
      engine_out.freq_buffer[n] =
        (beamformer[PYRAMIC_CHANNELS_IN * n + ch] * spectrum_in[PYRAMIC_CHANNELS_IN * n + ch]);
  }

  engine_out.synthesis(buffer_out);

  for (size_t n = 0 ; n < NUM_SAMPLES ; n++)
  {
    // adjust volume
    buffer_out[n] *= 8.;

    // clip
    if (buffer_out[n] > 1)
      output[2*n] = float2int;
    else if (buffer_out[n] < -1)
      output[2*n] = -float2int;
    else
      output[2*n] = (int16_t)(float2int * buffer_out[n]);

    // copy second channel
    output[2*n+1] = output[2*n];
  }

}

/*************************/
/*** THE INFINITE LOOP ***/
/*************************/

// Use this to exit the possibly infinite processing loop
bool is_running = true;
void signal_handler(int param)
{
  printf("The program was interrupted. Cleaning up now.");
  is_running = false;
}

// Now the main program
int main(void)
{
  int ret;
  float ellapsed_time = 0.;

  // install the signal handler
  std::signal(SIGINT, signal_handler);

  // User defined initializations
  init();

  // Start pyramic and the loop
  Pyramic pyramic(NUM_SAMPLES);

  ret = pyramic.start();

  if (ret)
  {
    while (is_running)  // The loop runs until catching a SIGINT (i.e. ctrl-C)
    {
      // wait for some samples to be available
      while (!pyramic.read_available())
        usleep(50);
      buffer_t &input_samples = pyramic.read_pop();

      // We get an output buffer
      if (!pyramic.play_available())
      {
        printf("Buffer underflow at processing");
        continue;
      }
      buffer_t &output_samples = pyramic.play_pop();

      // Call the processing routine
      processing(input_samples, output_samples);

      // Push back the buffers in the queues
      pyramic.read_push(input_samples);
      pyramic.play_push(output_samples);

      // move time forward
      ellapsed_time += (float)pyramic.n_samples() / pyramic.samplerate();
    }

    pyramic.stop();
  }
  else
  {
    printf("Failed to start Pyramic.");
  }
}
