%% Lyapunov Exponents with the Rosenstein's Algorithm
%% This code accompanies the powerpoint slides on Lyapunov exponents.

% Description
% - This activity will step through some basic processing steps for
%   calculating Lyapunov Exponents using Rosenstein's algorithm

%% 1) Lyapunov of a basic sinusoid
clear
clc

% Here we'll calculate the LyE of a sinusoid with exponential growth. We're
% adding the exponential growth to make sure the LyEs are positive and so
% the AMI function doesn't produce the artifacts it does with perfectly
% sinusoidal data.

% Sampling frequency
Fs = 1000;

% Time will serve as the independent variable
t = (0:1/Fs:40)';

% Frequency of the sinusoid
f = 1;

% The equation of a sinusoid.
y = exp(0.01*t).*sin(2*pi*f*t)

% plot the results
figure
plot(t,y)
xlabel('Time (s)')
ylabel('Y')

% Determining time lag and embedding dimension
% First we need to find the time lag. From experience we know the first
% minimum of a highly periodic function will be less than the period, so we
% use that many data point as an estimate for the maximum lag.
L = 1000;
[tau_s, ami_s] = AMI_Stergiou(y, L);

% Plot AMI function
figure
plot(ami_s(:,1),ami_s(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI';'Stergiou'})
grid minor

% Second we need the embedding dimension. A circle exists in 2 dimensions
% so the result here should not be surprising.
MaxDim = 12;
tau = 250; % this value is found from the AMI function above.

Rtol = 15;
Atol = 2;
speed = 1;
[dE, dim] = FNN(y, tau, MaxDim, Rtol, Atol, speed);

% Create a plot so we can see the results.
figure
plot(1:MaxDim, dim, 'o-k')
hold
plot(dE, dim(dE), 'or') % point where NN < 0.1 percent - overly restrictive?
xlabel('Embedding Dimension')
ylabel('% FNN')
grid on

% Plot state space
y_emb = embed(y,2,250)'

figure
plot(y_emb(:,1),y_emb(:,2),'.-')
xlabel('y(t)')
ylabel('y(t+\tau)')

% For Rosenstein's we need different parameters. The slope here works for
% the exponential sinusoid but may not for other time series. With this 
% theoretical data we know the period is exactly 1 s (f = 1 Hz and 1/f = 1 s) and that there are 
% always f*Fs data point per cycle. This and the MeanPeriod are in units of 
% seconds. Again you should get a small positive number.
tau = 250; % delay
m = 2; % embedding dim
slope = [0, 1, 10, 20];
MeanPeriod = 1; % in units of  seconds.
% Note that the next function takes about 37 s on my laptop to run. :)
tic
[LyE_RS, LyE_RL, out_R] = LyE_R(y, Fs, tau, m, slope, MeanPeriod, 1);
toc

%% 2) Application of Rosenstein algorithm to human data
clc
clear

% load ankle plantar and dorsiflexion data sampled at 60Hz
load('A4_LyE_hip_data.mat')
Fs = 60;

% plot the data
figure
plot(t,data)
xlabel('Time (s)')
ylabel('Hip angle (deg)')

% Determining time lag and embedding dimension
% First we need to find the time lag. From experience we know the first
% minimum of a highly periodic function will be less than the period, so we
% use that many data point as an estimate for the maximum lag.
L = 1000;
[tau_s, ami_s] = AMI_Stergiou(data, L);

% Plot AMI function
figure
plot(ami_s(:,1),ami_s(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI';'Stergiou'})
grid minor

% Second we need the embedding dimension. A circle exists in 2 dimensions
% so the result here should not be surprising.
MaxDim = 12;
tau = 20; % this value is found from the AMI function above.

Rtol = 15;
Atol = 2;
speed = 1;
[dE, dim] = FNN(data, tau(1), MaxDim, Rtol, Atol, speed);

% Create a plot so we can see the results.
figure
plot(1:MaxDim, dim, 'o-k')
hold
plot(dE, dim(dE), 'or') % point where NN < 0.1 percent - overly restrictive?
xlabel('Embedding Dimension')
ylabel('% FNN')
grid on

% Need to get an idea of the period of the attractor to set the evolve parameter:
% Visualize period of an attractor in 3D
data_emb = embed(data,4,20)';

figure
plot3(data_emb(1:200,1),data_emb(1:200,2),data_emb(1:200,3),'o-')
hold on
grid on
axis square

for smpl = 1:100
plot3(data_emb(1:1+smpl,1),data_emb(1:1+smpl,2),data_emb(1:1+smpl,3),'o-r','LineWidth',3)
title(['evole = ' num2str(smpl) ' You can ctrl+c to exit.'])
pause
end

% Another way to estaimate periodicity is to run a power spectruum estimate
[Pxx,F] = periodogram(data-mean(data),hamming(length(data)),length(data),Fs);

% Find the dominant frequency and then get the period as 1/Dom_freq and
% express in samples by multiplying by sampling rate (Fs)
period_fft_est = round((1/0.9111),2)

figure
plot(F,Pxx)
xlabel('Hz')
ylabel('Power (original units^2/Hz)')
title('Modified Periodogram Power Spectral Density Estimate')

% For Rosenstein's we need different parameters. Again you should get a
% small positive number.
m = 4; % embedding dimension
tau = 20; % time delay

slope = [0, .55, 10, 80]; % for gait data, a typical value for the first region is 0 to 1/2 period of attractor 
MeanPeriod = period_fft_est; % estimate of the period of the attractor (in seconds)
[LyE_RS, LyE_RL, out_R] = LyE_R(data, Fs, tau, m, slope, MeanPeriod, 1);
