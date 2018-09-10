import json, argparse
import numpy as np
import pyroomacoustics as pra

parser = argparse.ArgumentParser(description='Computes fixed beamforming weights for pyramic')
parser.add_argument('fs', type=int, help='The sampling frequency')
parser.add_argument('nfft', type=int, help='The length of the FFT')
parser.add_argument('mics_loc', type=str, help='Path to JSON file containing microphone locations')
parser.add_argument('output', type=str, help='The path to the output file')
args = parser.parse_args()

with open('pyramic.json', 'r') as f:
    data = json.load(f)
pyramic = np.array(data['pyramic'])
pyramic -= np.mean(pyramic, axis=1, keepdims=True)

t_room = 27  # degrees (temperature)
h_room = 58  # percent (humidity)
c = pra.parameters.calculate_speed_of_sound(t_room, h_room, 1000)

# a vector pointing from the source towards the array (horizontal plane)
p =  - (pyramic[7] - pyramic[0]) * np.r_[1,1,0]

delays = np.dot(pyramic, p) / c

omega = 2. * np.pi * (np.arange(nfft // 2 + 1) / nfft * fs)

# the weights are stored with channels interleaved: [ch0, ch1, ..., ch47, ...]
weights = np.exp(1j * omega[:,None] * delays[None,:]).flatten()

# store in complex interleaved format: [real(x0), imag(x0), real(x1), image(x1), ...]
data_out = {
        'fixed_weights' : np.c_[np.real(weights), np.imag(weights)].flatten().tolist(),
    }

with open(args.output, 'w') as f:
    json.dump(data_out, f)
