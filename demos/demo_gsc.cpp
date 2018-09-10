#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>

#include <gsc.h>
#include <stft.h>
#include <pyramic.h>

#define NUM_SAMPLES 1024
#define NFFT (2 * NUM_SAMPLES)

#define GSC_FILE_CONFIG "config/demo_gsc.json"
//#define GSC_FILE_WEIGHTS "config/demo_gsc_weights.json"
#define GSC_FILE_WEIGHTS "config/weights_delay_sum.json"

/*******************/
/* GLOBAL VARS ETC */
/*******************/

STFT engine_in(NUM_SAMPLES, NFFT, 0, 0, PYRAMIC_CHANNELS_IN, STFT_WINDOW_BOTH);
STFT engine_out(NUM_SAMPLES, NFFT, 0, 0, 1, STFT_WINDOW_BOTH);  // single output channel
GSC gsc(GSC_FILE_CONFIG, GSC_FILE_WEIGHTS, NFFT, PYRAMIC_SAMPLERATE, PYRAMIC_CHANNELS_IN);
float buffer_in[NUM_SAMPLES * PYRAMIC_CHANNELS_IN] = {0};
float buffer_out[NUM_SAMPLES] = {0};
float int2float = 1. / (1 << 15);
int16_t float2int = (1 << 15) - 1;

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init()
{
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
  // Convert the input buffer to float
  for (size_t n = 0 ; n < NUM_SAMPLES * PYRAMIC_CHANNELS_IN ; n++)
    buffer_in[n] = int2float * input[n];

  // Meat of the processing
  e3e_complex *spectrum_in = engine_in.analysis(buffer_in);
  gsc.process(spectrum_in, engine_out.freq_buffer);
  engine_out.synthesis(buffer_out);

  // Post-processing before output
  for (size_t n = 0 ; n < NUM_SAMPLES ; n++)
  {
    // clip and convert to int16_t
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
