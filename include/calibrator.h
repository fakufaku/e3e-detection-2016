#ifndef __CALIBRATOR_H__
#define __CALIBRATOR_H__

#include <string>

#include <Eigen/Dense>
#include <e3e_detection.h>

class Calibrator
{
  public:

    std::string output_filename;
    int nfft;
    int nchannel;
    int nfreq;
    int sample_count;
    Eigen::MatrixXcf covmat;   // Stack of covariance matrices (nchannel, nfreq * nchannel)
    
    Calibrator(std::string _output_filename, int _nfft, int _nchannel);
    void process(e3e_complex *input);
    void finalize();
    void save(std::string filename);
};

#endif // __CALIBRATOR_H__
