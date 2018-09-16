
#include <unistd.h>
#include <assert.h>
#include <rls.h>

RLS::RLS(int _nchannel, int _nfreq, float _ff, float _reg)
 : nchannel(_nchannel), nfreq(_nfreq), ff(_ff), reg(_reg)
{
  // for convenience
  this->ff_inv = 1.f / ff;
  
  // Initialize the buffers in the queue
  for (int n = 0 ; n < QSIZE ; n++)
  {
    this->input[n].setZero(this->nchannel, this->nfreq);
    this->ref_signal[n].setZero(this->nfreq);
    q_free.push(n);
  }

  // Initialize the weights
  this->weights[0].setZero(this->nchannel, this->nfreq);
  this->weights[1].setZero(this->nchannel, this->nfreq);

  // RLS variables
  this->covmat_inv.resize(this->nchannel, this->nfreq * this->nchannel);
  for (int f = 0 ; f < this->nfreq ; f++)
  {
    auto Rinv = this->covmat_inv.block(0, f * this->nchannel, this->nchannel, this->nchannel);
    Rinv.setIdentity(this->nchannel, this->nchannel);
    Rinv *= (1.f / this->reg);
  }
  this->xcov = Eigen::ArrayXXcf::Zero(this->nchannel, this->nfreq);

  // Now start the update thread
  this->is_running = true;
  this->run_thread = std::thread(&RLS::run, this);
}

RLS::~RLS()
{
  // stop the thread and wait until it terminates
  this->is_running = false;
  this->run_thread.join();
}

void RLS::push(Eigen::ArrayXXcf &new_input, Eigen::ArrayXcf &new_ref)
{
  int b = -1;  // the buffer index

  // If there is an unused buffer, use it
  this->mutex_q_free.lock();
  if (this->q_free.size() > 0)
  {
    b = this->q_free.front();
    this->q_free.pop();
  }
  this->mutex_q_free.unlock();

  // If we couldn't get an unused buffer, just use the oldest in the ready queue
  if (b < 0)
  {
    this->mutex_q_ready.lock();
    assert(this->q_ready.size() > 0);
    b = this->q_ready.front();
    this->q_ready.pop();
    this->mutex_q_ready.unlock();
  }

  // Fill the buffer with the data received
  this->input[b] = new_input;
  this->ref_signal[b] = new_ref;

  // add the buffer to ready queue
  this->mutex_q_ready.lock();
  this->q_ready.push(b);
  this->mutex_q_ready.unlock();
}

void RLS::fill(Eigen::ArrayXXcf &weights)
{
  this->mutex_weights.lock();
  weights = this->weights[this->w_consume_index];
  this->mutex_weights.unlock();
}

void RLS::update(Eigen::ArrayXXcf &input, Eigen::ArrayXcf &ref_signal)
{
  /*
   * Updates the inverse covariance matrix and cross-covariance vector.
   * Then, solves for the new adaptive weights
   *
   * @param input The input reference signal vector
   * @param error The error signal
   */
   float one_m_ff = (1. - this->ff);
   float ff_loc = this->ff / one_m_ff;

  // Update cross-covariance vector
  this->xcov = one_m_ff * (ff_loc * this->xcov + input.rowwise() * ref_signal.transpose().conjugate());

  // The rest needs to be done frequency wise
  for (int f = 0 ; f < this->nfreq ; f++)
  {
    // Update covariance matrix using Sherman-Morrison Identity
    auto Rinv = this->covmat_inv.block(0, f * this->nchannel, this->nchannel, this->nchannel);
    Eigen::VectorXcf u = Rinv * input.matrix().col(f);
    float v = 1. / (ff_loc + (input.matrix().col(f).adjoint() * u).real()(0,0)); // the denominator is a real number
    Rinv = this->ff_inv * (Rinv - (v * u * u.adjoint()));
    
    // Multiply the two to obtain the new adaptive weight vector
    this->weights[this->w_update_index].col(f) = Rinv * this->xcov.matrix().col(f);
  }
}

void RLS::run()
{
  int b;  // the buffer index

  while (this->is_running)
  {
    // pop a buffer from the ready for update list, if available
    this->mutex_q_ready.lock();
    if (this->q_ready.size() == 0)
    {
      // nothing to do, wait a bit and try again
      this->mutex_q_ready.unlock();
      usleep(200);
      continue;
    }
    else
    {
      b = this->q_ready.front();
      this->q_ready.pop();
      this->mutex_q_ready.unlock();
    }

    // compute the update
    this->update(this->input[b], this->ref_signal[b]);

    // put back the buffer in the queue
    this->mutex_q_free.lock();
    this->q_free.push(b);
    this->mutex_q_free.unlock();

    // swap the weights buffer
    this->mutex_weights.lock();
    int t = w_consume_index;
    w_consume_index = w_update_index;
    w_update_index = t;
    this->mutex_weights.unlock();
  }
}

