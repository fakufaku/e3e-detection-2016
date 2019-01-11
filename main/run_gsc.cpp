#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>

#include <AudioFile.h>
#include <json.hpp>
#include <hpf.h>
#include <gsc.h>
#include <stft.h>

using nlohmann::json;

/*******************/
/* GLOBAL VARS ETC */
/*******************/

STFT *engine_in, *engine_out;
GSC *gsc;
HPF *hpf;
float *buffer_in, *buffer_out;
int frame_size = 0;

// for the input audio file
AudioFile<float> afile_in;
int sample_index = 0;
int num_samples = 0;
int samplerate = 0;
int n_channels = 1;
int bitdepth = 0;

// for the input audio file
AudioFile<float> afile_out;
int sample_index_out;
int num_samples_out;
std::string afile_out_name = "temp.wav";

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init(int argc, char **argv)
{
  if (argc != 5)
  {
    printf("Usage: %s <config_file> <weights_file> <input> <output>\n", argv[0]);
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
  float rls_ff_inv = 1.f / rls_ff;              // ... and its inverse
  float rls_reg = config.at("rls_reg").get<float>();  // regularization factor for RLS
  float pb_ff = config.at("pb_ff").get<float>();      // forgetting factor for projection back
  int pb_ref_channel = config.at("pb_ref_channel").get<int>();  // The reference channel for projection back
  float high_pass_cutoff = config.at("high_pass_cutoff").get<float>();
  float f_max = config.at("f_max").get<float>();

  // Read in the input file
  afile_in.load(std::string(argv[3]));
  num_samples = afile_in.getNumSamplesPerChannel();
  samplerate = afile_in.getSampleRate();
  bitdepth = afile_in.getBitDepth();
  std::cout << "The file has " << num_samples << std::endl;

  // Prepare the output file
  afile_out_name = argv[4];
  // compute the number of output samples
  num_samples_out = (num_samples / frame_size) * frame_size;
  sample_index_out = 0;
  // Set both the number of channels and number of samples per channel
  afile_out.setAudioBufferSize (1, num_samples_out);
  afile_out.setBitDepth(bitdepth);
  afile_out.setSampleRate(samplerate);

  // allocate the sample buffers
  buffer_in = new float[frame_size * n_channels];
  buffer_out = new float[frame_size];

  // allocate the high pass filter
  hpf = new HPF(high_pass_cutoff, samplerate, nfft);

  // allocate GSC object
  gsc = new GSC(weights_file, nfft, samplerate, n_channels,
      nchannel_ds, rls_ff, rls_reg, pb_ff, pb_ref_channel, f_max);
  engine_in = new STFT(frame_size, nfft, 0, 0, n_channels, STFT_WINDOW_BOTH);
  engine_out = new STFT(frame_size, nfft, 0, 0, 1, STFT_WINDOW_BOTH);

  // Zero the DC and middle element of the output frequency buffer
  engine_out->freq_buffer[0] = 0.;
  engine_out->freq_buffer[nfft / 2] = 0.;
}

void clean_up()
{
  // save the output to file
  afile_out.save(afile_out_name, AudioFileFormat::Wave);

  delete gsc;
  delete hpf;
  delete engine_in;
  delete engine_out;
  delete[] buffer_in;
  delete[] buffer_out;
}

/************************************/
/* Read fresh samples from the file */
/************************************/

bool get_fresh_samples(float *buffer)
{
  /*
   * Fills the buffer with fresh samples. If no more are availables, 
   * it returns 0, otherwise 1.
   */

  for (int n = 0 ; n < frame_size ; n++, sample_index++)
  {
    if (sample_index >= num_samples)
    {
      for (int ch = 0 ; ch < n_channels ; ch++)
        buffer[n * n_channels + ch] = 0.;
    }
    else
    {
      for (int ch = 0 ; ch < n_channels ; ch++)
        buffer[n * n_channels + ch] = afile_in.samples[ch][sample_index];
    }
  }

  if (sample_index < num_samples)
    return true;
  else
    return false;
}

void push_new_samples(float *buffer)
{
  /*
   * Fills the output file buffer
   */

  for (int n = 0 ; n < frame_size && sample_index_out < num_samples_out ; n++)
    afile_out.samples[0][sample_index_out] = buffer[n];
}

/***********************************/
/* USER-DEFINED PROCESSING ROUTINE */
/***********************************/

// All the processing can be concentrated in this function
void processing(float *input, float *output)
{ 
  // These constants are needed to convert from int16 to float and back
  static float int2float = 1. / (1 << 15);
  static int16_t float2int = (1 << 15) - 1;

  // Convert the input buffer to float
  bool has_clipped_in = false;
  for (size_t n = 0 ; n < frame_size * n_channels ; n++)
  {
    if (!has_clipped_in && abs(input[n]) >= float2int)
    {
      std::cout << "Input: clipping" << std::endl;
      has_clipped_in = true;
    }
    buffer_in[n] = int2float * input[n];
  }

  // Go to freq domain
  e3e_complex *spectrum_in = engine_in->analysis(buffer_in);

  // generalized sidelobe canceller
  gsc->process(spectrum_in, engine_out->freq_buffer);

  // Apply HPF in frequency domain
  hpf->process(engine_out->freq_buffer);

  // go back to time domain
  engine_out->synthesis(buffer_out);

  // Post-processing before output
  bool has_clipped_out = false;
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
    output[2*n+1] = input[n_channels * n];

    // check for clipping
    if (!has_clipped_out && (abs(output[2*n]) >= float2int || abs(output[2*n+1] >= float2int)))
    {
      std::cout << "Output: clipping" << std::endl;
      has_clipped_out = true;
    }
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
  float ellapsed_time = 0.;

  // install the signal handler
  std::signal(SIGINT, signal_handler);

  // User defined initializations
  init(argc, argv);

  while (get_fresh_samples(buffer_in))  // The loop runs until catching a SIGINT (i.e. ctrl-C)
  {
    // Call the processing routine
    processing(buffer_in, buffer_out);

    // Push new samples to the output file
    push_new_samples(buffer_out);
  }

  clean_up();
}
