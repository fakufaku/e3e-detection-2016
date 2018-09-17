
#include <unistd.h>
#include <assert.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <rls.h>

RLS::RLS(int _nchannel, int _nfreq, float _ff, float _reg)
 : nchannel(_nchannel), nfreq(_nfreq), ff(_ff), reg(_reg)
{
  // for convenience
  this->ff_inv = 1.f / this->ff;
  this->one_m_ff = 1.f - this->ff;
  this->ff_ratio = this->ff / this->one_m_ff;
  this->reg_inv = 1.f / this->reg;
  
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
    std::cout << "RLS: Dropping frame" << std::endl;
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
  static int count = 0;
  count++;
  /*
   * Updates the inverse covariance matrix and cross-covariance vector.
   * Then, solves for the new adaptive weights
   *
   * @param input The input reference signal vector
   * @param error The error signal
   */
  // Update cross-covariance vector
  this->xcov = this->ff * this->xcov
             + this->one_m_ff * (input.rowwise() * ref_signal.transpose().conjugate());

  // The rest needs to be done frequency wise
  for (int f = 0 ; f < this->nfreq ; f++)
  {
    // Update covariance matrix using Sherman-Morrison Identity
    auto input_double_mat = input.cast<std::complex<double>>().matrix();
    auto Rinv = this->covmat_inv.block(0, f * this->nchannel, this->nchannel, this->nchannel);
    Eigen::VectorXcd u = Rinv * input_double_mat.col(f);

    // the denominator is a real number
    float v = 1. / ((double)this->ff_ratio + (input_double_mat.col(f).adjoint() * u).real()(0,0)); 

    // We know the diagonal elements should be purely real
    auto uuH = (v * u) * u.adjoint();

    Rinv = (double)this->ff_inv * (Rinv - ((v * u) * u.adjoint()));
    Rinv.diagonal() = Rinv.diagonal().real();  // we now this should be real and positive

    if (f == 50 && count == 100)
    {
      std::cout << "The trace is " << Rinv.trace() << std::endl;
      count = 0;
    }
    
    // Multiply the two to obtain the new adaptive weight vector
    this->weights[this->w_update_index].col(f) = (Rinv * this->xcov.cast<std::complex<double>>().matrix().col(f)).array().cast<std::complex<float>>();
  }

  // swap the weights buffer
  this->mutex_weights.lock();
  int t = this->w_consume_index;
  this->w_consume_index = this->w_update_index;
  this->w_update_index = t;
  this->mutex_weights.unlock();
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
  }
}

