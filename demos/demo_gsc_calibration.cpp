#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cmath>
#include <complex>
#include <iostream>

#include <e3e_detection.h>
#include <pyramic.h>
#include <AudioFile.h>
#include <stft.h>
#include <calibrator.h>

/*******************/
/* GLOBAL VARS ETC */
/*******************/

Calibrator *calib;
float *buffer_in, *buffer_out;
int frame_size = 0;
int nfft = 0;
float calibration_duration = 0.;

float buffer_in[NUM_SAMPLES * PYRAMIC_CHANNELS_IN] = {0};
float buffer_out[NUM_SAMPLES] = {0};
float int2float = 1. / (1 << 15);
int16_t float2int = (1 << 15) - 1;

// We will first store all the frames, then process them
std::vector<e3e_complex_vector> all_frames;

/*****************************/
/* USER-DEFINED INIT ROUTINE */
/*****************************/

void init(int argc, char **argv)
{
  if (argc != 4)
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
  nfft = config.at("nfft").get<int>();

  // allocate the sample buffers
  buffer_in = new float[frame_size * PYRAMIC_CHANNELS_IN];
  buffer_out = new float[frame_size];

  engine_in = new STFT(frame_size, nfft, 0, 0, PYRAMIC_CHANNELS_IN, STFT_WINDOW_BOTH);

  // initialize the calibration module
  calib = new Calibrator(weights_file, nfft, PYRAMIC_CHANNELS_IN);
}

void clean_up()
{
  // Now process all the collected frames
  std::cout << "Start processing the data collected" << std::endl;
  for (int i = 0 ; i < all_frames.size() ; i++)
    calib.process(all_frames[i].data());
  calib.finalize();
  std::cout << "Processing finished, cleaning up" << std::endl;

  delete buffer_in;
  delete buffer_out;
  delete calib;
  delete engine_in;
}

/***********************************/
/* USER-DEFINED PROCESSING ROUTINE */
/***********************************/

// All the processing can be concentrated in this function
void processing(std::vector<float> &input, buffer_t &output)
{ 
  // These constants are needed to convert from int16 to float and back
  static float int2float = 1. / (1 << 15);
  static int16_t float2int = (1 << 15) - 1;

  // Convert the input buffer to float
  for (size_t n = 0 ; n < frame_size * PYRAMIC_CHANNELS_IN ; n++)
    buffer_in[n] = int2float * input[n];

  // Spectral analysis
  e3e_complex *spectrum_in = engine_in->analysis(buffer_in);

  // now copy the spectrum to a new std vector
  // This is not very efficient, but we probably don't care all that much...
  int frame_total_size = (nfft / 2 + 1) * PYRAMIC_CHANNELS_IN;
  std::vector<e3e_complex_vector> arr(frame_total_size);
  for (int i = 0 ; i < frame_total_size ; i++)
    arr[i] = spectrum_in[i];
  all_frames.push_back(arr);

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
