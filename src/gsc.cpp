#include <fstream>
#include <cmath>
#include <assert.h>
#include <iostream>

#include <json.hpp>
#include <gsc.h>

using nlohmann::json;

GSC::GSC(
      std::string weights_file,           // path to the file containing the fixed beamforming weights
      int _nfft,     // the FFT size
      float _fs,      // the sampling frequency
      int _nchannel,   // the number of input channels
      int _nchannel_ds,
      float _rls_ff,
      float _rls_reg,
      float _pb_ff,
      int _pb_ref_channel,
      float _f_max
    )
    : nfft(_nfft), fs(_fs), nchannel(_nchannel), nchannel_ds(_nchannel_ds),
      rls_ff(_rls_ff), rls_reg(_rls_reg), pb_ff(_pb_ff),
      pb_ref_channel(_pb_ref_channel), f_max(_f_max)
{
  /**
    Constructor of the Generalized Sidelobe Canceller (GSC) class.

    @param nfft The FFT size used on the input signal
    @param fs The sampling frequency of the input signal
    @param nchannel The number of input channels
    @param nchannel_ds The number of channels after the spatial downsampling step
    @param rls_ff The forgetting factor of RLS
    @param rls_reg The regularization parameter of RLS
    @param pb_ff The forgetting factor of the projection back
    @param pb_ref_channel The reference input channel for the projection back
    @param f_max The maximum frequency to process, above this frequency everything is set to zero
    */

  // compute the downsampling factor
  this->ds = this->nchannel / this->nchannel_ds;
  this->ds_inv = 1.f / this->ds;

  // Limit frequencies
  this->f_min_index = 1;  // we skip the DC component in the processing
  this->f_max_index = int(ceilf((this->f_max / this->fs) * this->nfft + 0.5)); // round to closest bin
  this->nfreq = this->f_max_index - this->f_min_index;  // only consider the number of bands processed

  // Read the file that contains the weights
  std::ifstream f_weights(weights_file);
  json j_weights;
  f_weights >> j_weights;
  f_weights.close();

  std::cout << "Finished reading weights" << std::endl << std::flush;

  // Get the fixed weights from the json file, the complex numbers are stored
  // with real/imag parts interleaved i.e. [r0, i0, r1, i1, r2,  ...]
  // in row-major order
  std::vector<float> w = j_weights.at("fixed_weights").get<std::vector<float>>();
  assert((int)w.size() == 2 * (this->nfft / 2 + 1) * this->nchannel);  // check the size is correct
  this->fixed_weights.setZero(this->nchannel, this->nfreq);
  for (int f = 0, offset = this->f_min_index * this->nchannel ; f < this->nfreq ; f++, offset += this->nchannel)
  {
    for (int ch = 0 ; ch < this->nchannel ; ch++)
      this->fixed_weights(ch, f) = e3e_complex(w[2 * (offset + ch)], w[2 * (offset + ch) + 1]);
  }
  // Normalize the weights for the projection
  this->fixed_weights.matrix().colwise().normalize();
  
  // Size the other buffers as needed
  this->adaptive_weights.setZero(this->nchannel_ds, this->nfreq);

  // Intermediate buffers
  this->output_fixed.setZero(this->nfreq);
  this->output_blocking_matrix.setZero(this->nchannel, this->nfreq);
  this->input_adaptive.setZero(this->nchannel_ds, this->nfreq);

  // Projection back buffers
  this->projback_num.setOnes(this->nfreq);
  this->projback_den.setOnes(this->nfreq);

  // RLS variables
  this->covmat_inv.resize(this->nchannel_ds, this->nfreq * this->nchannel_ds);
  for (int f = 0 ; f < this->nfreq ; f++)
  {
    auto Rinv = this->covmat_inv.block(0, f * this->nchannel_ds, this->nchannel_ds, this->nchannel_ds);
    Rinv.setIdentity(this->nchannel_ds, this->nchannel_ds);
    Rinv *= (1.f / this->rls_reg);
  }

  this->xcov = Eigen::ArrayXXcf::Zero(this->nchannel_ds, this->nfreq);

  std::cout << "Finished initialization of GSC object." << std::endl << std::flush;
}

void GSC::process(e3e_complex *input, e3e_complex *output)
{
  // Pre-emptivaly zero-out the content of output buffer
  for (int f = 0 ; f < this->nfft / 2 + 1 ; f++)
    output[f] = 0.f;

  // Wrap input/output in Eigen::Array
  // WARNING: Eigen library stores arrays/matrices column-major (by default)
  int input_offset = this->f_min_index * this->nchannel;
  Eigen::Map<Eigen::ArrayXXcf> X(input + input_offset, this->nchannel, this->nfreq);
  Eigen::Map<Eigen::ArrayXcf> Y(output + this->f_min_index, this->nfreq);

  // Compute the fixed beamformer output
  this->output_fixed = (this->fixed_weights.conjugate() * X).colwise().sum();

  // Apply the blocking matrix
  this->output_blocking_matrix = X - this->fixed_weights.rowwise() * this->output_fixed.transpose();

  // Downsample the channels to a reasonnable number
  for (int c = 0, offset = 0 ; c < this->nchannel_ds ; c++, offset += this->ds)
    this->input_adaptive.row(c) = this->output_blocking_matrix.block(offset, 0, this->ds, this->nfreq).colwise().sum() * this->ds_inv;

  // Update the adaptive weights
  this->rls_update(input_adaptive, this->output_fixed);

  // Compute the output signal
  Y = this->output_fixed - (this->adaptive_weights.conjugate() * this->input_adaptive).colwise().sum().transpose();

  // projection back: apply scale to match the output to channel 1
  this->projback(X, Y, this->pb_ref_channel);

}

void GSC::rls_update(Eigen::ArrayXXcf &input, Eigen::ArrayXcf &ref_signal)
{
  /*
   * Updates the inverse covariance matrix and cross-covariance vector.
   * Then, solves for the new adaptive weights
   *
   * @param input The input reference signal vector
   * @param error The error signal
   */
   float one_m_ff = (1. - this->rls_ff);
   float ff = this->rls_ff / one_m_ff;

  // Update cross-covariance vector
  this->xcov = one_m_ff * (ff * this->xcov + input.rowwise() * ref_signal.transpose().conjugate());

  // The rest needs to be done frequency wise
  for (int f = 0 ; f < this->nfreq ; f++)
  {
    // Update covariance matrix using Sherman-Morrison Identity
    auto Rinv = this->covmat_inv.block(0, f * this->nchannel_ds, this->nchannel_ds, this->nchannel_ds);
    Eigen::VectorXcf u = Rinv * input.matrix().col(f);
    float v = 1. / (ff + (input.matrix().col(f).adjoint() * u).real()(0,0)); // the denominator is a real number
    Rinv = this->rls_ff_inv * (Rinv - (v * u * u.adjoint()));
    
    // Multiply the two to obtain the new adaptive weight vector
    this->adaptive_weights.col(f) = Rinv * this->xcov.matrix().col(f);
  }
}

void GSC::projback(Eigen::Map<Eigen::ArrayXXcf> &input, Eigen::Map<Eigen::ArrayXcf> &output, int input_ref_channel)
{
  /*
   * This function updates the projection back weight and scales
   * the output with the new coefficient
   */

  // slice out the chosen columns of input
  this->projback_num = this->pb_ff * this->projback_num + (1.f - this->pb_ff) * (output.conjugate() * input.row(input_ref_channel).transpose());
  this->projback_den = this->pb_ff * this->projback_den + (1.f - this->pb_ff) * output.abs2();

  // weight the output
  output *= (this->projback_num / this->projback_den);
}

