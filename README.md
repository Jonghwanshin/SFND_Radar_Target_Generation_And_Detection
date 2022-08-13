# Radar Target Generation and Detection

## Project Overview

This project aims to determine range and velocity of an target object from given radar signal. 

This includes:

- Configure the FMCW waveform based on the system requirements.
- Define the range and velocity of target and simulate its displacement.
- For the same simulation loop process the transmit and receive signal to determine the beat signal
- Perform Range FFT on the received signal to determine the Range
- Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.

![project overview](images/project-overview.png)

## Radar System Requirements

- In this project, I designed a Radar based on the given system requirements (above).

![Radar System Requirements](images/radar-req.png)



## Writeup

### FP01. FMCW Waveform Design

I could determine the parameters (B, Tchirp, slope) for designing a FMCW waveform by below equations
```math
B_{sweep} = c / 2 * res\\
T_{chirp} = 5.5 * 2 * R_{max} / c\\
Slope = B / T_{chirp}
```
The result of slope is `2.0455e+13` which meets the specifications.

### FP02. Simulation Loop

I modeled the transmitted signal and received signal with  below equations:
```math
Tx(t) = cos(2\pi(f_c t + \frac{\alpha * t ^2}{2})) \\
Rx(t) = cos(2\pi(f_c (t-\tau) + \frac{\alpha * (t-\tau) ^2}{2}))
```
The time delay $\tau$ is the travel time of transmitted signal therefore I modeled $\tau$ as below equations:
```matg
\tau = 2 * R(t)/c \\
R(t) = R_{init} + V_{init} \times t;
```

### FP03. Range FFT (1st FFT)

I implemented 1D FFT for the given radar measurement with below MATLAB code.

```MATLAB
%% RANGE MEASUREMENT
Mix_1d = reshape(Mix, [Nr, Nd]);  % reshape the vector into Nr*Nd array. 
                                  % Nr and Nd defines the size of
                                  % Range and Doppler FFT respectively.
Mix_fft = fft(Mix_1d,Nr)/Nr; % run the FFT on the beat signal 
                             % along the range bins dimension (Nr) and
                             % normalize.
Mix_fft = abs(Mix_fft);      % Take the absolute value of FFT output
Mix_fft = Mix_fft(1:Nr/2); % throw out the half since it is double sided.
						   % it will also leave the first FFT result only.
```

I set the initial position of a target object to 110m and the velocity to -20m/s in this project and I can observe the target object in the below FFT results.

| 1D FFT                       | 2D FFT                       |
| ---------------------------- | ---------------------------- |
| ![1D FFT](images/1d-fft.png) | ![2D FFT](images/2d-fft.png) |

### FP04. 2D CFAR

I implemented 2D CFAR with `conv2`, <u>which is 2D convolution function</u> of MATLAB with below procedures.

1. decide the training cell(`Td, Tr`) and guard cell(`Gd, Gr`) offsets and threshold in dB(`threshold`).

2. create a convolution kernel to sum all values from training cells.
3. perform 2D convolution to Radar Doppler Map(RDM) that converted into power and divide into # of training cells to get noise level.
4. calculate `threshold` from noise level  and `offset` .
5. pad `threshold` and convert into dB scale and compare to RDM.
6. plot result.

```MATLAB
% Step1. decide parameters for 2D CFAR
Tr = 10; Td = 4; % Training Cells in both dimensions
Gr = 5; Gd = 2; % Guard Cells in both dimensions
offset = 3; % offset the threshold by SNR value in dB

% Step2. create conv kernel
kernel_cfar = ones(2*Tr+2*Gr+1, 2*Td+2*Gd+1);
cut_r = Tr+Gr+1; % middle of kernel
cut_d = Td+Gd+1;
kernel_cfar(cut_r-Gr:cut_r+Gr,cut_d-Gd:cut_d+Gd) = 0; % remove guard cells
size_kernel = sum(kernel_cfar, 'all'); % get number of training cells

% Step 3~5. get noise_level from RDM, calculate threshold and compare to RDM
threshold = offset*(conv2(db2pow_(RDM), kernel_cfar, "valid")/size_kernel);
pad_r = round(size(kernel_cfar,1)/2-1);
pad_d = round(size(kernel_cfar,2)/2-1);
threshold = padarray(pow2db_(threshold), [pad_r, pad_d], inf);
is_above_threshold = double(RDM > threshold);

% Step 6. plot result
figure,surf(doppler_axis,range_axis,is_above_threshold);
colorbar;

% define db2pow and pow2db 
% since I don't have Signal Processing Toolbox
function pow = db2pow_(db)
    pow = 10.^(db/10);
end

function db = pow2db_(pow)
    db = 10 * log10(pow);
end
```

The output of 2D CFAR is below graph and this describes the target object I defined.

![2D CFAR](images/2d-cfar.png)
