#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>

#include <pyramic.h>
#include <AudioFile.h>
#include <stft.h>
#include <calibrator.h>

#define NUM_SAMPLES 1024
#define NFFT (2 * NUM_SAMPLES)

/*******************/
/* GLOBAL VARS ETC */
/*******************/

float buffer_in[NUM_SAMPLES * PYRAMIC_CHANNELS_IN] = {0};
float buffer_out[NUM_SAMPLES] = {0};
float int2float = 1. / (1 << 15);
int16_t float2int = (1 << 15) - 1;

Calibrator calib("data/the_weights.json", NFFT, PYRAMIC_CHANNELS_IN);

// for the audio file
AudioFile<float> afile;
int sample_index = 0;
int num_samples = 0;

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init(int argc, char **argv)
{
  std::string fn(argv[1]);
  afile.load(fn);
  num_samples = afile.getNumSamplesPerChannel();
  std::cout << "The file has " << num_samples << std::endl;
}

void clean_up()
{
  calib.finalize();
}

/************************************/
/* Read fresh samples from the file */
/************************************/

bool get_fresh_samples(std::vector<float> &buffer)
{
  /*
   * Fills the buffer with fresh samples. If no more are availables, 
   * it returns 0, otherwise 1.
   */

  for (int n = 0 ; n < NUM_SAMPLES ; n++, sample_index++)
  {
    if (sample_index >= num_samples)
    {
      for (int ch = 0 ; ch < PYRAMIC_CHANNELS_IN ; ch++)
        buffer[n * PYRAMIC_CHANNELS_IN + ch] = 0.;
    }
    else
    {
      for (int ch = 0 ; ch < PYRAMIC_CHANNELS_IN ; ch++)
        buffer[n * PYRAMIC_CHANNELS_IN + ch] = afile.samples[ch][sample_index];
    }
  }

  if (sample_index < num_samples)
    return true;
  else
    return false;
}

/***********************************/
/* USER-DEFINED PROCESSING ROUTINE */
/***********************************/

// All the processing can be concentrated in this function
void processing(std::vector<float> &input, buffer_t &output)
{ 
  static STFT engine_in(NUM_SAMPLES, NFFT, 0, 0, PYRAMIC_CHANNELS_IN, STFT_WINDOW_BOTH);

  // Convert the input buffer to float
  for (size_t n = 0 ; n < NUM_SAMPLES * PYRAMIC_CHANNELS_IN ; n++)
    buffer_in[n] = input[n];

  // Meat of the processing
  e3e_complex *spectrum_in = engine_in.analysis(buffer_in);
  calib.process(spectrum_in);

  // Post-processing before output
  for (size_t n = 0 ; n < NUM_SAMPLES ; n++)
  {
    // copy second channel
    output[2*n] = 0;
    output[2*n+1] = 0;
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
int main(int argc, char **argv)
{
  int ret;
  std::vector<float> input_samples(NUM_SAMPLES * PYRAMIC_CHANNELS_IN);

  // install the signal handler
  std::signal(SIGINT, signal_handler);

  // User defined initializations
  init(argc, argv);

  // Start pyramic and the loop
  Pyramic pyramic(NUM_SAMPLES);

  ret = pyramic.start();

  if (ret)
  {
    while (is_running)  // The loop runs until catching a SIGINT (i.e. ctrl-C)
    {
      // wait for some samples to be available
      is_running = get_fresh_samples(input_samples);

      if (!is_running)
        break;

      // Wait for free output buffer
      while (is_running && !pyramic.play_available())
        usleep(50);

      buffer_t &output_samples = pyramic.play_pop();

      // Call the processing routine
      processing(input_samples, output_samples);

      // Push back the buffers in the queues
      pyramic.play_push(output_samples);
    }

    pyramic.stop();
  }
  else
  {
    printf("Failed to start Pyramic.");
  }

  clean_up();
}
