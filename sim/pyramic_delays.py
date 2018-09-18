import json, argparse
import numpy as np
import matplotlib.pyplot as plt
import pyroomacoustics as pra

type_choices = ['ds', 'ds-mean']

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Computes fixed beamforming weights for pyramic')
    parser.add_argument('fs', type=int, help='The sampling frequency')
    parser.add_argument('nfft', type=int, help='The length of the FFT')
    parser.add_argument('mics_loc', type=str, help='Path to JSON file containing microphone locations')
    parser.add_argument('output', type=str, help='The path to the output file')
    parser.add_argument('-t', '--type', type=str, choices=type_choices, default=type_choices[0],
            help='The type of beamformer to compute')
    parser.add_argument('-b', '--beta', type=float, default=60.,
            help='Target angle of the beamformer')
    args = parser.parse_args()

    with open(args.mics_loc, 'r') as f:
        data = json.load(f)

    pyramic = np.array(data['pyramic'])  # shape : (nmics, ndim)
    pyramic -= np.mean(pyramic, axis=0, keepdims=True)
    pyramic_pc = pra.experimental.PointCloud(X=pyramic.T)

    t_room = 27  # degrees (temperature)
    h_room = 58  # percent (humidity)
    c = pra.parameters.calculate_speed_of_sound(t_room, h_room, 1000)

    # Let us consider a source that is about 70cm away with center a few centimeters below the center of the array
    source_loc = np.r_[0.7, 0., -0.05,]

    # a vector pointing from the source towards the array (horizontal plane)
    p =  source_loc / np.linalg.norm(source_loc)

    # the frequency vector
    omega = 2. * np.pi * (np.arange(args.nfft // 2 + 1) / args.nfft * args.fs)

    delays = np.dot(pyramic, p) / c
    ds_weights = np.exp(-1j * omega[:,None] * delays[None,:])  # shape: (nfreq, nmics)

    if args.type == 'ds':
        weights = ds_weights

    elif args.type == 'ds-mean':

        g = pra.doa.GridSphere(n_points=300)
        beta_half = np.radians(args.beta) / 2
        ip = 1. - np.dot(g.cartesian.T, p)
        nz = ip <= beta_half
        all_delays = np.dot(pyramic, g.cartesian[:,nz]) / c  # shape: (nmics, n_steering_vecs)
        steering_vectors = np.exp(-1j * omega[:,None,None] * all_delays[None,:,:])  # shape: (nfreq, nmics, n_steering_vecs)
        steering_vectors /= np.linalg.norm(steering_vectors, axis=1, keepdims=True)

        print(steering_vectors.shape[2])

        weights = np.zeros(steering_vectors.shape[:2], dtype=steering_vectors.dtype)

        for f, omg in enumerate(omega):
            u,s,v = np.linalg.svd(steering_vectors[f,:,:])
            weights[f,:] = u[:,0]


        # plot the target response
        g.values = nz.astype(np.int)
        g.plot(plotly=False)
        plt.title('Desired response')

    
    # Plot a few response
    gp = pra.doa.GridSphere(n_points=1000)
    all_delays = np.dot(pyramic, gp.cartesian) / c  # shape: (nmics, n_steering_vecs)
    for f_hz in [500., 1000., 2000., 4000., 8000., 12000, 16000]:

        f = np.argmin(np.abs(omega - 2 * np.pi * f_hz))
        steering_vectors = np.exp(-1j * 2* np.pi * f_hz * all_delays)  # shape: (nmics, n_steering_vecs)
        steering_vectors /= np.linalg.norm(steering_vectors, axis=1, keepdims=True)
        gp.values = 20. * np.log10(np.abs(np.dot(steering_vectors.T, np.conj(weights[f,:]))))
        gp.values = np.abs(np.dot(steering_vectors.T, np.conj(weights[f,:]))) ** 2
        gp.plot(plotly=False)
        plt.title('{} Hz'.format(f_hz))

    plt.show()


    # the weights are stored with channels interleaved: [ch0, ch1, ..., ch47, ...]
    weights_flat = weights.flatten()

    # store in complex interleaved format: [real(x0), imag(x0), real(x1), image(x1), ...]
    data_out = {
            'fixed_weights' : np.c_[np.real(weights_flat), np.imag(weights_flat)].flatten().tolist(),
        }

    with open(args.output, 'w') as f:
        json.dump(data_out, f)
