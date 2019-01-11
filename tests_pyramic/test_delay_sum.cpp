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
 1.55827705e-04,  8.80705102e-05,  2.03133149e-05, -6.78956315e-06,
-2.03410022e-05, -4.74438803e-05, -1.15201076e-04, -1.82958271e-04,
-2.09217221e-04, -1.62467186e-04, -1.15717151e-04, -9.70171365e-05,
-8.76671294e-05, -6.89671154e-05, -2.22170802e-05,  2.45329550e-05,
 1.91130177e-04,  1.51614960e-04,  1.12099742e-04,  9.62936547e-05,
 8.83906112e-05,  7.25845241e-05,  3.30693063e-05, -6.44591141e-06,
-9.11424055e-05, -1.19384383e-04, -1.47626360e-04, -1.58923151e-04,
-1.64571547e-04, -1.75868338e-04, -2.04110315e-04, -2.32352293e-04,
 2.14265249e-04,  1.93258089e-04,  1.72250929e-04,  1.63848065e-04,
 1.59646633e-04,  1.51243769e-04,  1.30236609e-04,  1.09229449e-04,
 5.98354268e-05,  4.13273692e-05,  2.28193115e-05,  1.54160884e-05,
 1.17144769e-05,  4.31125384e-06, -1.41968038e-05, -3.27048615e-05
};

/*************************/
/* CREATE THE BEAMFORMER */
/*************************/

void compute_beamforming_weights()
{
  const std::complex<float> j(0, 1);  // imaginary number
  float denom = 1. / PYRAMIC_CHANNELS_IN;

  for (size_t f = 0 ; f < (NFFT / 2 + 1) ; f++)
  {
    double f_hz = (double)f / NFFT * PYRAMIC_SAMPLERATE;
    for (size_t ch = 0 ; ch < PYRAMIC_CHANNELS_IN ; ch++)
      beamformer[f * PYRAMIC_CHANNELS_IN + ch] = std::exp(-j * (float)(2. * M_PI * f_hz * delays[ch])) * denom;
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
    // now apply delay and sum the other channels
    for (size_t ch = 0 ; ch < PYRAMIC_CHANNELS_IN ; ch++)
      engine_out.freq_buffer[n] =
        (std::conj(beamformer[PYRAMIC_CHANNELS_IN * n + ch]) * spectrum_in[PYRAMIC_CHANNELS_IN * n + ch]);
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
