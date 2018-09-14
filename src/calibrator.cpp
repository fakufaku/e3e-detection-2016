#include <fstream>
#include <iostream>
#include <json.hpp>
#include <calibrator.h>

// for convenience
using json = nlohmann::json;

Calibrator::Calibrator(std::string _output_filename, int _nfft, int _nchannel)
 : output_filename(_output_filename), nfft(_nfft), nchannel(_nchannel)
{
  this->nfreq = this->nfft / 2 + 1;
  this->sample_count = 0;
  covmat.setZero(this->nfreq * this->nchannel, this->nchannel);
}

void Calibrator::process(e3e_complex *input)
{
  Eigen::Map<Eigen::MatrixXcf> X(input, this->nfreq, this->nchannel);
  this->sample_count++;

  for (int f = 1 ; f < this->nfreq ; f++)  // skip DC (start at 1)
  {
    auto R = this->covmat.block(f * this->nchannel, 0, this->nchannel, this->nchannel);
    auto x = X.row(f).transpose();
    R += x * x.adjoint();
  }
}

void Calibrator::finalize()
{
  // linear vector of weights with interleaved real/imag parts
  std::vector<float> weights; // will be of size (2 * this->nfreq * this->nchannel) at the end

  // fill in the missing DC component
  for (int n = 0 ; n < 2 * this->nchannel ; n++)
    weights.push_back(0.);

  // now process the rest
  for (int f = 1 ; f < this->nfreq ; f++)  // skip DC (start at 1)
  {
    auto R = this->covmat.block(f * this->nchannel, 0, this->nchannel, this->nchannel);
    R /= this->sample_count;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcf> eigensolver(R);
    if (eigensolver.info() != Eigen::Success)
    {
      // Print an error message and fill the weights with junk
      std::cout << "Error: the EVD did not converge for f = " << f << std::endl;
      for (int n = 0 ; n < 2 * this->nchannel ; n++)
        weights.push_back(0.);
      continue;
    }

    // find the location of largest eigenvalue
    int max_ev;
    eigensolver.eigenvalues().maxCoeff(&max_ev);

    // and the corresponding eigenvector
    auto w = eigensolver.eigenvectors().col(max_ev);

    // Add in the weight vector in interleaved fashion
    for (int n = 0 ; n < this->nchannel ; n++)
    {
      weights.push_back(std::real(w(2 * n)));
      weights.push_back(std::imag(w(2 * n + 1)));
    }
  }

  // Finally, open a json file and write the weights in it
  json j;
  j["fixed_weights"] = weights;
  std::ofstream o(this->output_filename);
  o << j;
  o.close();
}
