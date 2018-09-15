#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>

#include <gsc.h>
#include <stft.h>
#include <pyramic.h>

/*******************/
/* GLOBAL VARS ETC */
/*******************/

STFT *engine_in, *engine_out;
GSC *gsc;
float *buffer_in, *buffer_out;
int frame_size = 0;

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init(int argc, char **argv)
{
  if (argc != 3)
  {
    printf("Usage: %s <config_file> <weights_file>\n", argv[0]);
    exit(0);
  }

  // The configuration files
  std::string config_file(argv[1]);
  std::string weights_file(argv[2]);
  
  // read in the JSON file containing the configuration
  std::ifstream i(config_file, std::ifstream::in);
  json config;
  i >> config;
  i.close();

  std::cout << "Finished reading config file" << std::endl << std::flush;

  // get all the GSC parameters
  frame_size = config.at("frame_size").get<int>();
  int nfft = config.at("nfft").get<int>();
  int nchannel_ds = config.at("nchannel_ds").get<int>();
  float rls_ff = config.at("rls_ff").get<float>();    // forgetting factor for RLS
  float rls_ff_inv = 1.f / this->rls_ff;              // ... and its inverse
  float rls_reg = config.at("rls_reg").get<float>();  // regularization factor for RLS
  float pb_ff = config.at("pb_ff").get<float>();      // forgetting factor for projection back
  int pb_ref_channel = config.at("pb_ref_channel").get<int>();  // The reference channel for projection back
  float f_max = config.at("f_max").get<float>();

  // allocate the sample buffers
  buffer_in = new float[frame_size * PYRAMIC_CHANNELS_IN];
  buffer_out = new float[frame_size];

  // allocate GSC object
  gsc = new GSC(weights_file, nfft, PYRAMIC_SAMPLERATE, PYRAMIC_CHANNELS_IN,
      nchannels_ds, rls_ff, rls_reg, pb_ff, pb_ref_channel, f_max);
  engine_in = new STFT(frame_size, nfft, 0, 0, PYRAMIC_CHANNELS_IN, STFT_WINDOW_BOTH);
  engine_out = new STFT(frame_size, nfft, 0, 0, 1, STFT_WINDOW_BOTH);

  // Zero the DC and middle element of the output frequency buffer
  engine_out->freq_buffer[0] = 0.;
  engine_out->freq_buffer[NFFT / 2] = 0.;
}

void clean_up()
{
  delete gsc;
  delete engine_in;
  delete engine_out;
  delete buffer_in;
  delete buffer_out;
}

/***********************************/
/* USER-DEFINED PROCESSING ROUTINE */
/***********************************/

// All the processing can be concentrated in this function
void processing(buffer_t &input, buffer_t &output)
{ 
  // These constants are needed to convert from int16 to float and back
  static float int2float = 1. / (1 << 15);
  static int16_t float2int = (1 << 15) - 1;

  // Convert the input buffer to float
  for (size_t n = 0 ; n < frame_size * PYRAMIC_CHANNELS_IN ; n++)
    buffer_in[n] = int2float * input[n];

  // Meat of the processing
  e3e_complex *spectrum_in = engine_in->analysis(buffer_in);
  gsc->process(spectrum_in, engine_out->freq_buffer);
  engine_out->synthesis(buffer_out);

  // Post-processing before output
  for (size_t n = 0 ; n < frame_size ; n++)
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
  Pyramic pyramic(frame_size);

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

  clean_up();
}
