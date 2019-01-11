#include "pyramic.h"

uint32_t toggle_half(uint32_t half)
{ 
  if (half == 1)
    return 2;
  else
    return 1;
}


Pyramic::Pyramic(int _num_samples)
 : num_samples(_num_samples)
{
  int i;

  this->read_buffer_size = PYRAMIC_CHANNELS_IN * this->num_samples;
  this->play_buffer_size = PYRAMIC_CHANNELS_OUT * this->num_samples;

  for (i = 0 ; i < N_BUFFERS ; i++)
  {
    this->read_buffers[i].resize(this->read_buffer_size);
    this->q_read_empty.push(&(this->read_buffers[i]));
    this->play_buffers[i].resize(this->play_buffer_size);
    this->q_play_empty.push(&(this->play_buffers[i]));
  }
}

Pyramic::~Pyramic()
{
  // clean up
  if (this->flag_is_running || this->p != NULL)
    this->stop();
}

int Pyramic::start()
{
  // Initialize the Pyramic array
  this->p = pyramicInitializePyramic();
  
  if (p != NULL)
  {

    this->flag_is_running = true;

    // start reader and player threads
    this->play_thread = std::thread(&Pyramic::player, this);
    this->read_thread = std::thread(&Pyramic::reader, this);


    return 1;  // success
  }
  else
  {
    return 0;  // failure
  }
}

void Pyramic::stop()
{
  this->flag_is_running = false;

  // Wait for threads to terminate
  this->read_thread.join();
  this->play_thread.join();

  // Now turn off pyramic
  if (this->p != NULL)
    pyramicDeinitPyramic(this->p);
}

void Pyramic::reader()
{
  uint32_t i;
  uint32_t current_half = 1;
  int16_t *input_buffers[2];

  // Start capture
  pyramicStartCapture(this->p, 2 * this->num_samples);
  struct inputBuffer inBuf = pyramicGetInputBuffer(this->p, 0); // 0 for 1st half, actually a pointer to the whole buffer

  // the double buffer
  input_buffers[0] = inBuf.samples;
  input_buffers[1] = inBuf.samples + this->read_buffer_size;

  // start the infinite loop
  while (this->flag_is_running)
  {
    // Wait while current buffer is busy
    while(pyramicGetCurrentBufferHalf(p) == current_half)
    {
      if (!this->flag_is_running)
        goto stop;
      usleep(50);
    }

    this->mutex_read_empty.lock();
    if (!this->q_read_empty.empty())
    {
      // get an empty buffer from the queue
      buffer_t &buffer = *(this->q_read_empty.front());
      this->q_read_empty.pop();
      this->mutex_read_empty.unlock();

      // fill it
      for (i = 0 ; i < this->read_buffer_size ; i++)
        buffer[i] = input_buffers[current_half-1][i];

      // put it in the other queue
      this->mutex_read_ready.lock();
      this->q_read_ready.push(&buffer);
      this->mutex_read_ready.unlock();
    }
    else
    {
      this->mutex_read_empty.unlock();
      printf("Buffer overflow at reading\n");
    }

    current_half = toggle_half(current_half);
  }

stop:

  pyramicStopCapture(this->p);
}

void Pyramic::player()
{

  uint32_t i;
  uint32_t current_half = 1;

  // make sure the playback module is off at the beginning
  pyramicEnableOutput(this->p, 0);

  // Get pointer to output buffer
  struct outputBuffer outBuf = pyramicGetOutputBuffer(this->p, 4 * this->play_buffer_size);  // 2nd arg is length in bytes (i.e., 2 buffers x 2 bytes/sample)

  // The double buffer
  int16_t *out_buffers[2] = {
    outBuf.samples,
    outBuf.samples + this->play_buffer_size
  };

  // zero out the first buffer
  for (i = 0 ; i < outBuf.length ; i++)
    outBuf.samples[i] = 0;

  // configure the playback module
  pyramicSelectOutputSource(this->p, SRC_MEMORY);  // use samples from memory
  pyramicSetOutputBuffer(this->p, outBuf);  // configure the buffer pointer and size

  // wait for a couple of buffers to accumulate
  while (this->q_play_ready.size() < PLAY_MIN_BUF_IN_Q)
    usleep(50);

  // start the loop
  while (this->flag_is_running)
  {
    // Wait for current half to be idle
    while (pyramicGetCurrentOutputBufferHalf(this->p) == current_half)
    {
      if (!this->flag_is_running)
        goto stop;
      usleep(50);
    }

    this->mutex_play_ready.lock();
    if (!this->q_play_ready.empty())
    {
      // get a buffer ready to play
      buffer_t &buffer = *(this->q_play_ready.front());
      this->q_play_ready.pop();
      this->mutex_play_ready.unlock();

      // copy it to the output buffer
      for (i = 0 ; i < this->play_buffer_size ; i++)
        out_buffers[current_half-1][i] = buffer[i];

      // put back the used buffer in the queue
      this->mutex_play_empty.lock();
      this->q_play_empty.push(&buffer);
      this->mutex_play_empty.unlock();

      // enable the Pyramic playback module if needed
      if (pyramicGetCurrentOutputBufferHalf(this->p) == 0)
      {
        pyramicEnableOutput(this->p, 1);

        // wait for the module to start
        while (pyramicGetCurrentOutputBufferHalf(this->p) == 0)
          usleep(50);
      }
    }
    else
    {
      this->mutex_play_ready.unlock();
      // fill with zeros when nothing is available
      for (i = 0 ; i < this->play_buffer_size ; i++)
        out_buffers[current_half-1][i] = 0;
      printf("Buffer underflow at playback\n");
    }

    current_half = toggle_half(current_half);
  }

stop:
    
  // zero out the output buffer at the end
  for (i = 0 ; i < outBuf.length ; i++)
    outBuf.samples[i] = 0;

  // disable the playback module
  pyramicEnableOutput(this->p, 0);

}
