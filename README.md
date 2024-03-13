# Antenna Array Beamforming — MATLAB implementation


### What's beamforming?
In a conventional beamformer, each antenna‘s signal is multiplied by a specific Fourier weighting, and the weighted signals are then summed together. 


<img src="beamforming scheme.png" width="700"/>

The weights are chosen to optimize the array’s response in a particular direction, forming a “beam” of heightened sensitivity in that direction. 

This process allows **precise detection and localization** of transmitted and received signals.


### My project is divided into 4 parts

1. Array mathematics
2. Direction of arrival (DOA) estimation
3. QPSK modulation and demodulation
4. Simulation result and application
