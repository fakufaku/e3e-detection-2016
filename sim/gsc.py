import numpy as np
from scipy import linalg
import pyroomacoustics as pra

from utilities import Plottable, ShortTimeAverage, LeakyIntegration

def calibration(signal):
    '''
    This routine will compute the fixed beamforming weights
    from the STFT of a recorded calibration signal.

    Parameters
    ----------
    signal: array_like (nframes, nfrequencies, nchannels)
        The STFT of the recorded calibration signal

    Returns
    -------
    array_like (nfrequencies, nchannels)
        The beamforming weights (in frequency domain)
    '''

    nframes, nfreq, nchan = signal.shape

    weights = np.zeros(signal.shape[1:], dtype=signal.dtype)

    for f in range(nfreq):

        R = np.dot(signal[:,f,:].T, np.conj(signal[:,f,:])) / nframes
        val, vec = linalg.eigh(R, eigvals=(R.shape[0]-1,R.shape[0]-1))
        weights[f,:] = vec[:,0]

    # normalize
    weights /= np.linalg.norm(weights, axis=1, keepdims=True)

    return weights


class GSC(Plottable):

    def __init__(self, calibration_signal, step_size, pb_ff, nfft, shift=None, win=None):

        self.step_size = step_size  # for the LMS
        self.pb_ff = pb_ff  # forgetting factor for the projection back
        self.nfft = nfft
        self.nchannel = calibration_signal.shape[1]

        if shift is None:
            shift = nfft // 2

        if win is None:
            win = pra.hann(nfft)

        # Compute the fixed beamformer
        X = pra.transform.analysis(calibration_signal, nfft, shift, win=win)
        self.fixed_weights = calibration(X)[1:,:]  # remove DC

        # gsc adaptive weights
        self.adaptive_weights = np.zeros_like(self.fixed_weights)

        # projection back weights
        self.pb_den = np.ones(self.fixed_weights.shape[0], dtype=self.fixed_weights.dtype)
        self.pb_num = np.ones(self.fixed_weights.shape[0], dtype=np.float)

    def process(self, X):

        assert X.shape[0] == self.nfft // 2 + 1
        assert X.shape[1] == self.nchannel
        output = np.zeros(self.nfft // 2 + 1, dtype=X.dtype)

        # ignore DC
        X = X[1:,:]

        # the fixed branch of the GSC
        out_fixed_bf = np.sum(np.conj(self.fixed_weights[1:,:]) * X, axis=1)

        # Now run the online projection back (seems more natural to do it on the fixed branch)
        self.pb_den = ( self.pb_ff * self.pb_den 
                + (1 - self.pb_ff) * np.conj(out_fixed_bf) * X[:,0] )
        self.pb_num = ( self.pb_ff * self.pb_num
                + (1 - self.pb_ff) * np.conj(out_fixed_bf) * out_fixed_bf )
        out_fixed_bf *= (self.pb_den / self.pb_num)

        # the adaptive branch of the GSC
        noise_ss_signal = X - self.fixed_weights * out_fixed_bf[:,None]  # projection onto null space
        noise_ss_norm = np.sum(noise_ss_signal * np.conj(noise_ss_signal), axis=1)
        out_adaptive_bf = np.sum(np.conj(self.adaptive_weights) * noise_ss_signal, axis=1)

        # compute the error signal and update the beamformer
        output[1:] = out_fixed_bf - out_adaptive_bf
        mu = (self.step_size / noise_ss_norm[:,None])
        self.adaptive_weights +=  mu * noise_ss_signal * output[1:,None]

        return output

class GSC_Newton(Plottable):

    def __init__(self, calibration_signal, step_size, pb_ff, nfft, ds, shift=None, win=None):

        self.step_size = step_size  # for the LMS
        self.pb_ff = pb_ff  # forgetting factor for the projection back
        self.nfft = nfft
        self.nchannel = calibration_signal.shape[1]
        self.ds = ds

        if shift is None:
            shift = nfft // 2

        if win is None:
            win = pra.hann(nfft)

        # Compute the fixed beamformer
        X = pra.transform.analysis(calibration_signal, nfft, shift, win=win)
        self.fixed_weights = calibration(X)[1:,:]  # remove DC
        self.norm_weights = 1. / np.linalg.norm(self.fixed_weights, axis=1, keepdims=True) ** 2

        # gsc adaptive weights
        self.adaptive_weights = np.zeros(
                (self.fixed_weights.shape[0], self.nchannel // ds),
                dtype=self.fixed_weights.dtype
                )
        self.adaptive_weights[:,0] = 1.

        # projection back weights
        self.pb_den = np.ones(self.fixed_weights.shape[0], dtype=self.fixed_weights.dtype)
        self.pb_num = np.ones(self.fixed_weights.shape[0], dtype=np.float)

        self.estimates = {
                'covmat' : LeakyIntegration(
                    0.9,  # average over this number of frames
                    lambda X : X[:,:,None] * np.conj(X[:,None,:]),  # (nfreq, nchan, nchan),
                    init=np.array([np.eye(self.nchannel // self.ds) for i in range(self.nfft // 2)]) * 1e-3,
                    ),
                'xcov' : LeakyIntegration(
                    0.9,
                    lambda v : v[0] * np.conj(v[1][:,None]),
                    ),
                }


    def process(self, X):

        assert X.shape[0] == self.nfft // 2 + 1
        assert X.shape[1] == self.nchannel
        output = np.zeros(self.nfft // 2 + 1, dtype=X.dtype)

        # remove DC
        X = X[1:,:]

        # the fixed branch of the GSC
        out_fixed_bf = np.sum(np.conj(self.fixed_weights) * X, axis=1)

        # the adaptive branch of the GSC
        noise_ss_signal = (X - self.norm_weights * self.fixed_weights * out_fixed_bf[:,None]).reshape((-1,self.nchannel // self.ds, self.ds)).sum(axis=-1)  # projection onto null space, and downsampling of channels

        # update covariance matrix estimate
        self.estimates['covmat'].update(noise_ss_signal)
        self.estimates['xcov'].update([noise_ss_signal, out_fixed_bf])
        covmat = self.estimates['covmat'].get()
        xcov = self.estimates['xcov'].get()

        ad_w = np.array([np.linalg.solve(covmat[n] + 1e-5, xcov[n]) 
            for n in range(covmat.shape[0])])

        out_adaptive_bf = np.sum(np.conj(ad_w) * noise_ss_signal, axis=1)

        # compute the error signal and update the beamformer
        output[1:] = out_fixed_bf - out_adaptive_bf

        # Now run the online projection back (seems more natural to do it on the fixed branch)
        self.pb_num = ( self.pb_ff * self.pb_num 
                + (1 - self.pb_ff) * np.conj(output[1:]) * X[:,0] )
        self.pb_den = ( self.pb_ff * self.pb_den
                + (1 - self.pb_ff) * np.conj(output[1:]) * output[1:] )
        output[1:] *= (self.pb_num / self.pb_den)

        #mu = self.step_size / (np.sum(whitened * np.conj(whitened), axis=1) + 1e-10)
        #self.adaptive_weights =  0.9 * self.adaptive_weights  + 0.1 * ad_w
        self.adaptive_weights = ad_w

        return output
