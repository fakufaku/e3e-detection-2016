#ifndef __GSC_H__
#define __GSC_H__
/*
 * This is the class definition for the Generalized Sidelobe Canceller
 * 
 * 2018 (c) Robin Scheibler
 */

#include <string>

#include <Eigen/Dense>
#include <e3e_detection.h>

class GSC
{
  public:


    // parameter attributes
    int nfft;
    float fs;
    int nchannel;

    // parameters coming from config file
    int nchannel_ds;  // number of channels after downsampling in adaptive branch
    int ds;           // number of channels combined together for downsampling (ds = nchannel / nchannel_ds)
    float ds_inv;     // == 1. / ds
    int nfreq;
    
    // algorithm parameters
    float rls_ff;
    float rls_ff_inv;
    float rls_reg;
    float pb_ff;
    int pb_ref_channel;

    // The limit of the processing band in frequency
    float f_max;
    int f_min_index, f_max_index;

    // The beamforming weights
    Eigen::ArrayXXcf fixed_weights;    // size: (nfreq, nchannel)
    Eigen::ArrayXXcf adaptive_weights; // size: (nfreq, nchannel_ds)

    // The intermediate buffers
    Eigen::ArrayXcf output_fixed;             // size: (nfreq)
    Eigen::ArrayXXcf output_blocking_matrix;  // size: (nfreq, nchannels)
    Eigen::ArrayXXcf input_adaptive;          // size: (nfreq, nchannels_ds)

    // Projection back variables
    Eigen::ArrayXcf projback_num;  // numerator, size: (nfreq), complex
    Eigen::ArrayXf projback_den;  // denominator, size: (nfreq), real-valued

    // RLS variables
    std::vector<Eigen::MatrixXcf> covmat_inv;  // inverse covariance matrices, size: (nfreq, nchannel_ds, nchannel_ds)
    Eigen::ArrayXXcf xcov;                     // cross covariance vectors, size: (nfreq, nchannel_ds)

    GSC(
        std::string fixed_beamformer_file,  // path to the configuration file
        std::string weights_file,           // path to the file containing the fixed beamforming weights
        int _nfft,     // the FFT size
        float fs,      // the sampling frequency
        int nchannel   // the number of input channels
       );

    void process(e3e_complex_vector &input, e3e_complex_vector &output);    
    void rls_update(Eigen::Map<Eigen::ArrayXXcf> &input, Eigen::ArrayXcf &ref_signal);
    void projback(Eigen::Map<Eigen::ArrayXXcf> &input, Eigen::Map<Eigen::ArrayXcf> &output, int input_ref_channel);
};

#endif // __GSC_H__
