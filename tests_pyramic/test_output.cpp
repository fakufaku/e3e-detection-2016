/*
 * main.c
 *
 *  Created on: Dec 22, 2016
 *      Author: ferry
 */


#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <thread>
#include <queue>
#include <cmath>

#include "pyramicio.h"

#define BUFFER_LEN 100
#define SINE_PERIOD 48

int sine_counter = 0;
int16_t sine[SINE_PERIOD];

// The Pyramic stuff
struct pyramic *p = NULL;

// Create a queue to store output frames
int is_playing = false;

uint32_t toggle_half_output(uint32_t half)
{ 
  if (half == 1)
    return 2;
  else
    return 1;
}

void playback()
{
  uint32_t i;
  int ret = 0;
  bool is_pyramic_enabled = false;

  // Get pointer to output buffer
  pyramicEnableOutput(p, 0);
  struct outputBuffer outBuf = pyramicGetOutputBuffer(p, 8 * BUFFER_LEN);  // 2nd arg is length in bytes (i.e., 2 buffers x 2 bytes/sample x 2 channels x buffer length)
  int16_t *out_buffers[2] = {
    outBuf.samples,
    outBuf.samples + 2 * BUFFER_LEN  // 1 buffer x 1 int16_t/sample x 2 channel x buffer length
  };

  // zero out the first buffer
  for (i = 0 ; i < 2 * 2 * BUFFER_LEN ; i++)
    outBuf.samples[i] = 0;

  uint32_t current_half = 1;

  // mark as currently running
  is_playing = true;

  // configure the output
  pyramicSelectOutputSource(p, SRC_MEMORY);
  ret = pyramicSetOutputBuffer(p, outBuf);
  if (ret == 0)
  {

    printf("Pyramic status: %d\n", (int)pyramicGetCurrentOutputBufferHalf(p));

    // start the loop
    while (is_playing)
    {
      // Wait for current half to be idle
      while (is_pyramic_enabled && pyramicGetCurrentOutputBufferHalf(p) == current_half)
        usleep(50);

      //printf("Pyramic current half: %d\n", (int)pyramicGetCurrentOutputBufferHalf(p));

      for (i = 0 ; i < BUFFER_LEN ; i++)
      {
        out_buffers[current_half-1][2*i] = sine[sine_counter];
        out_buffers[current_half-1][2*i+1] = sine[sine_counter];

        // increment sine_counter
        if (sine_counter == SINE_PERIOD - 1)
          sine_counter = 0;
        else
          sine_counter++;
      }

      if (not is_pyramic_enabled)
      {
        pyramicEnableOutput(p, 1);
        is_pyramic_enabled = true;
        // wait for pyramic to start
        while (pyramicGetCurrentOutputBufferHalf(p) == 0)
          usleep(50);
      }

      current_half = toggle_half_output(current_half);
    }
    
    // Stop output
    pyramicEnableOutput(p, 0);

    // zero out the output buffer at the end
    for (i = 0 ; i < 2 * 2 * BUFFER_LEN ; i++)
      outBuf.samples[i] = 0;
  }

}


int main(void)
{
  uint32_t i;
  double Ts = 2. * M_PI / SINE_PERIOD;

  // fill the sine wave
  for (i = 0 ; i < SINE_PERIOD ; i++)
  {
    sine[i] = (int16_t)(0.02 * (2 << 15) * sin(i * Ts));
  }

  // Initialize the Pyramic array
  p = pyramicInitializePyramic();

	if(p) {
		printf("Success in Initializing Pyramic !\n");

    // Start the playback too
    std::thread playback_thread(playback);

    sleep(20);

    is_playing = false;
    playback_thread.join();

	}
	else
		printf("Failed to init Pyramic !\n");

  pyramicDeinitPyramic(p);

	return 0;
}
