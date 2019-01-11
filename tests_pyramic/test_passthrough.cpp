/*
 * main.c
 *
 *  Created on: Dec 22, 2016
 *      Author: ferry
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <thread>
#include <queue>

#include "pyramicio.h"

#define NBUFFERS 3
#define BUFFER_LEN 2048
#define OCHANNELS 2
#define OUT_BUFFER_SIZE (OCHANNELS * BUFFER_LEN)

// The Pyramic stuff
struct pyramic *p = NULL;

// Create a queue to store output frames
std::queue<int16_t *> q_ready;  // store buffers ready to be played
std::queue<int16_t *> q_empty;  // store unused buffers
int16_t buffers[NBUFFERS][OUT_BUFFER_SIZE];
int is_playing = false;

uint32_t toggle_half_pt(uint32_t half)
{ 
  if (half == 1)
    return 2;
  else
    return 1;
}

void playback()
{
  uint32_t i;
  bool is_pyramic_enabled = false;

  // Get pointer to output buffer
  struct outputBuffer outBuf = pyramicGetOutputBuffer(p, 8 * BUFFER_LEN);  // 2nd arg is length in bytes (i.e., 2 buffers x 2 bytes/sample x 2 channels x buffer length)
  int16_t *out_buffers[2] = {
    outBuf.samples,
    outBuf.samples + 2 * BUFFER_LEN  // 1 buffer x 1 int16_t/sample x 2 channel x buffer length
  };

  // zero out the first buffer
  for (i = 0 ; i < outBuf.length ; i++)
    outBuf.samples[i] = 0;

  int16_t *buffer;
  uint32_t current_half = 1;

  // mark as currently running
  is_playing = true;

  // configure the output
  pyramicSelectOutputSource(p, SRC_MEMORY);
  pyramicSetOutputBuffer(p, outBuf);

  // wait for a couple of buffers to accumulate
  while (q_ready.size() < NBUFFERS - 1)
    usleep(50);

  // start the loop
  while (is_playing)
  {
    // Wait for current half to be idle
    while (is_pyramic_enabled && pyramicGetCurrentOutputBufferHalf(p) == current_half)
      usleep(50);

    if (!q_ready.empty())
    {
      buffer = q_ready.front();
      q_ready.pop();
      for (i = 0 ; i < OUT_BUFFER_SIZE ; i++)
        out_buffers[current_half-1][i] = buffer[i];
      q_empty.push(buffer);

      if (not is_pyramic_enabled)
      {
        pyramicEnableOutput(p, true);
        is_pyramic_enabled = true;

        while (pyramicGetCurrentOutputBufferHalf(p) == 0)
          usleep(50);
      }
    }
    else
    {
      for (i = 0 ; i < OUT_BUFFER_SIZE ; i++)
        out_buffers[current_half-1][i] = 0;
      printf("Buffer underflow at playback\n");
    }

    current_half = toggle_half_pt(current_half);
  }
    
  // zero out the output buffer at the end
  for (i = 0 ; i < outBuf.length ; i++)
    outBuf.samples[i] = 0;

}


int main(void)
{
  uint32_t i;
  uint32_t current_half = 1;
  int16_t *input_buffers[2];
  int16_t *buffer;

  // Initialize the Pyramic array
  p = pyramicInitializePyramic();

	if(p) {
		printf("Success in Initializing Pyramic !\n");

    // Fill the queue of empty buffers
    for (i = 0 ; i < NBUFFERS ; i++)
      q_empty.push(buffers[i]);

    // Start capture
		pyramicStartCapture(p, 2 * BUFFER_LEN);
    struct inputBuffer inBuf = pyramicGetInputBuffer(p, 0); // 0 for 1st half, actually a pointer to the whole buffer
    input_buffers[0] = inBuf.samples;
    input_buffers[1] = inBuf.samples + inBuf.microphoneCount * BUFFER_LEN;

    // Start the playback too
    std::thread playback_thread(playback);

    int n = 0;
    int nmax = (int)(20. * 48000 / BUFFER_LEN);
    while (n < nmax)
    {
      // Wait while current buffer is busy
      while(pyramicGetCurrentBufferHalf(p) == current_half)
        usleep(50);

      if (!q_empty.empty())
      {
        buffer = q_empty.front();
        q_empty.pop();
        for (i = 0 ; i < BUFFER_LEN ; i++)
        {
          buffer[2*i] = input_buffers[current_half-1][48*i];
          buffer[2*i+1] = input_buffers[current_half-1][48*i+1];
        }
        q_ready.push(buffer);
      }
      else
      {
        printf("Buffer overflow at reading\n");
      }

      current_half = toggle_half_pt(current_half);

      n++;
    }

    is_playing = false;
    playback_thread.join();

		pyramicStopCapture(p);
	}
	else
		printf("Failed to init Pyramic !\n");

  pyramicDeinitPyramic(p);

	return 0;
}
