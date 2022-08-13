clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 77e9;
res = 1;
r_max = 200;
v_max = 100;
c = 3e8;
%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
% You will provide the initial range and velocity of the target. 
% Range cannot exceed the max value of 200m and velocity can be any value 
% in the range of -70 to + 70 m/s. 
R_init = 110;
V_init = -20;


%% FMCW Waveform Generation

% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of 
% the FMCW chirp using the requirements above.
B = c / (2 * res);
Tchirp = 5.5 * 2 * r_max / c;
slope = B/Tchirp;
                                                
% The number of chirps in one sequence. 
% Its ideal to have 2^ value for the ease of running the FFT
% for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods 
                          % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    % For each time stamp 
    % update the Range of the Target for constant velocity. 
    r_t(i) = R_init + V_init * t(i);
    td(i) = 2 * r_t(i)/c;

    % update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + slope*t(i)^2/2)); 
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + slope*(t(i)-td(i))^2/2));

    Mix(i) = Tx(i) * Rx(i); % mix Tx and Rx signals
end

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
%plotting the range
figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 
plot(Mix_fft);
 
axis ([0 200 0 1]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10; Td = 4;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5; Gd = 2;

% offset the threshold by SNR value in dB
offset = 3;

% it slides the CUT across range doppler map by giving margins 
% at the edges for Training and Guard Cells.
% For every iteration sum the signal level within all the training cells.
% To convert the value from logarithmic to linear using db2pow function. 
% Average the summed values for all of the training cells used. 
% After averaging convert it back to logarithimic using pow2db.
% Further add the offset to it to determine the threshold. 
% Next, compare the signal under CUT with this threshold. 
% If the CUT level > threshold assign it a value of 1, else equate it to 0.

% define conv kernel
kernel_cfar = ones(2*Tr+2*Gr+1, 2*Td+2*Gd+1);
cut_r = Tr+Gr+1; % middle of kernel
cut_d = Td+Gd+1;
kernel_cfar(cut_r-Gr:cut_r+Gr,cut_d-Gd:cut_d+Gd) = 0;
size_kernel = sum(kernel_cfar, 'all'); % get number of training cells

% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
threshold = offset*(conv2(db2pow_(RDM), kernel_cfar, "valid")/size_kernel);
pad_r = round(size(kernel_cfar,1)/2-1);
pad_d = round(size(kernel_cfar,2)/2-1);
threshold = padarray(pow2db_(threshold), [pad_r, pad_d], inf);
is_above_threshold = double(RDM > threshold);

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
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