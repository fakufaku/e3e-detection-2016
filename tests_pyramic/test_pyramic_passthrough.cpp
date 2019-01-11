#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include "pyramic.h"

#define NUM_SAMPLES 1024

/***************************/
/* CORE PROCESSING ROUTINE */
/***************************/

// All the processing can be concentrated in this function
void processing(buffer_t &input, buffer_t &output,
                size_t n_samples, size_t channels_in,
                size_t channels_out, size_t samplerate)
{
  // Copy the first two input channels to the output
  for (size_t n = 0 ; n < n_samples ; n++)
  {
    output[channels_out * n]     = input[channels_in * n];
    output[channels_out * n + 1] = input[channels_in * n + 1];
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
      processing(input_samples, output_samples,
                 pyramic.n_samples(), pyramic.channels_in(),
                 pyramic.channels_out(), pyramic.samplerate());

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
