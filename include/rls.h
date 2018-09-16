#ifndef __RLS_H__
#define __RLS_H__

#include <Eigen/Dense>
#include <queue>
#include <thread>
#include <mutex>

#define QSIZE 2

class RLS
{
  public:
    int nchannel;
    int nfreq;
    float ff;
    float reg;
    float ff_inv;

    bool is_running = false;  // indicates wether the thread has been started or not
    
    Eigen::ArrayXXcf input[QSIZE];
    Eigen::ArrayXcf ref_signal[QSIZE];
    std::queue<int> q_ready;  // store buffers ready to be played
    std::mutex mutex_q_ready;  // protects the queue
    std::queue<int> q_free;  // store buffers ready to be played
    std::mutex mutex_q_free;  // protects the queue

    Eigen::ArrayXXcf weights[2];  // use two buffers to quickly swap
    int w_consume_index = 0;
    int w_update_index = 1;
    std::mutex mutex_weights;  // protects the queue


    // the threads
    std::thread run_thread;

    // RLS variables
    Eigen::MatrixXcf covmat_inv;  // inverse covariance matrices, size: (nchannel, nnfreq * channel_ds), all matrices are stacked (horizontally)
    Eigen::ArrayXXcf xcov;        // cross covariance vectors, size: (nchannel, nfreq)

    RLS(int nchannel, int nfreq, float ff, float reg);
    ~RLS();

    void push(Eigen::ArrayXXcf &new_input, Eigen::ArrayXcf &new_ref);
    void fill(Eigen::ArrayXXcf &weights);
    void update(Eigen::ArrayXXcf &input, Eigen::ArrayXcf &ref_signal);
    void run();
};

#endif // __RLS_H__
