clear variables;
close all;
%% Load in signal vector
data = importdata('distanceForFFT.mat')

%% Calculate the FFT of the signal
signalFFT = fft(double(data));
P2 = abs(signalFFT/length(data));
P1 = P2(1:length(data)/2 +1);
P1(2:end-1) = 2*P1(2:end-1);

%% Define frequency domain and plot single sided amplitude spectrum p1
f = (1/20) * (0:(length(data)/2))/length(data);
plot(f,P1);