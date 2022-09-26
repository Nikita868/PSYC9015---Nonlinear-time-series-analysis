%% Lyapunov Exponents with Wulf Algorithm
%% This code accompanies the powerpoint slides on Lyapunov exponents.

% Description
% - This activity will step through some basic processing steps for
%   calculating Lyapunov Exponents using Wulf's algorithm

%% 1 Lyapunov of a basic sinusoid
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

% Wulf's Algorithm

% Input parameters:
tau = 250; % selected delay from AMI
m = 2; % selected embedding dimension from FNN
evolve = 0.3*Fs; % evolution between replacements (in samples). Fs is sampling rate (in Hz)

% Run Wulf's estimate:
[out_w, lyew] = LyE_W_NK(y, Fs, tau, m, evolve,0);

% Plot segment LyE estimate and the final estimate
figure
plot(out_w(:,6))
xlabel('Segment')
ylabel('LyE estimate')
hold
line([0 max(out_w(:,1))],[lyew lyew],'color','r')

lyew

%% 2 Choosing the "evolve" parameter

% What value of evolve should I choose?
% Run Wulf's estimate for various values of evolve parameter to see what
% you get. 

lyew_tot=[]
for evolve = 100:10:400
    
[out_w, lyew] = LyE_W(y, Fs, tau, m, evolve,0);
lyew_tot = [lyew_tot; lyew]

end

figure
plot(100:10:400,lyew_tot,'o-')
xlabel('Evolve parameter (samples)')
ylabel('LyE estimate')

%% 3) Compare simulated sine waves
% Simulate the same 1Hz sine wave with different rates of exponenetial
% growth, and use the Wulf algorithm to estimate the largest LyE of the
% time series. Use three values of the with exponential function  r1 =
% 0.01, r2 = 0.02, and r3 = 0.03. What do you think would happen to the
% estimate as you increase r? What do you notice about the LyE estimate
% from the Wulf algorithm?

y = exp(0.01*t).*sin(2*pi*f*t);
[out_w, lyew1] = LyE_W_NK(y, Fs, tau, m, evolve,0);
lyew1

y = exp(0.02*t).*sin(2*pi*f*t);
[out_w, lyew2] = LyE_W_NK(y, Fs, tau, m, evolve,0);
lyew2

y = exp(0.03*t).*sin(2*pi*f*t);
[out_w, lyew3] = LyE_W_NK(y, Fs, tau, m, evolve,1);
lyew3

%% 3.1 Potential problem in Wulf algorithm application:
% Make sure your reference and comparison trajectories are not on the same
% attractor cycle! This problem is illustrated in this simulation selection
% of a neighboring point that lies on the same trajectory.

% Let's see how this problem looks:

% Sampling frequency
Fs = 1000;

% Time will serve as the independent variable
t = (0:1/Fs:40)';

% Frequency of the sinusoid
f = 1;

% simulate new data
% The equation of a sinusoid.
y = exp(0.02*t).*sin(2*pi*f*t)

% plot the results
figure
plot(t,y)
xlabel('Time (s)')
ylabel('Y')

% Wulf's Algorithm

% Input parameters:
tau = 250; % selected delay from AMI
m = 2; % selected embedding dimension from FNN
evolve = 0.3*Fs; % evolution between replacements (in samples). Fs is sampling rate (in Hz)

% Run Wulf's estimate:
[out_w, lyew] = LyE_W_NK(y, Fs, tau, m, evolve,1);


%% 4) Application of Wulf algorithm to human data

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
title(['evole = ' num2str(smpl)])
pause
end

% Another way to estaimate periodicity is to run a power spectruum estimate
[Pxx,F] = periodogram(data-mean(data),hamming(length(data)),length(data),Fs);

figure
plot(F,Pxx)
xlabel('Hz')
ylabel('Power (original units^2/Hz)')
title('Modified Periodogram Power Spectral Density Estimate')

%  Wulf algorithm estimate
period_samples = 55; % based on 4x of first local min AMI value, but could also base on FFT result

tau = 20; % selected delay from AMI
m = 4; % selected embedding dimension from FNN
evolve = round(0.3*period_samples); % evolution between replacements (in samples). Fs is sampling rate (in Hz)

% Run Wulf's estimate for a small subset of the data to visualize the algorithm:
[out_w, lyew] = LyE_W_NK(data(1:200), Fs, tau, m, evolve,1);

% Run Wulf's estimte on the whole dataset
[out_w, lyew] = LyE_W_NK(data, Fs, tau, m, evolve,0);

% Plot segment LyE estimate and the final estimate
figure
plot(out_w(:,6))
xlabel('Segment')
ylabel('LyE estimate')
hold
line([0 max(out_w(:,1))],[lyew lyew],'color','r')

lyew

%% Choosing the evolve value
% In systems where the mechanism for chaos is unknown, one must check for
% exponent stability over a wide range of evolution try running over a
% range of evolve values, choose from a plateu

for j = 1:56
    [out_w, lyew] = LyE_W_NK(data, Fs, tau, m, j, 0);    
    lyew_tot(j) = lyew 
    j
end

figure
plot(lyew_tot,'o-')
xlabel('Evolve parameter (samples)')
ylabel('LyE estimate')

%%  5) Application of Wulf algorithm to the Lorenz attractor data
[t,y] = ChaosLibrary('Lorenz',[0,100],[0.1 -0.01 9],[10 28 8/3]);

% remove initial transient
y = y(500:end,:)
t = t(500:end,:)

% plot individual variables of Lorenz
figure
subplot(311)
plot(t,y(:,1))
ylabel('X')

subplot(312)
plot(t,y(:,2))
ylabel('Y')

subplot(313)
plot(t,y(:,3))
ylabel('Z')
xlabel('Time (s)')

% Average Mutual Information and Choosing time delay for embedding

% This will run the AMI_Stergiou method (histogram) with a max lag of 150
% on the first state variable y(:,1) of the Rossler system. The choice of
% maximum lag is largely determined through experience.
L = 150;
[tau_s, ami_s] = AMI_Stergiou(y(:,1),L);

% Plot AMI function
figure
plot(ami_s(:,1),ami_s(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI';'Stergiou'})
grid minor

% FNN
MaxDim = 12;
tau = 12; % this value is found from the AMI function above.

Rtol = 15;
Atol = 2;
speed = 0;
[dE,dim] = FNN(y(:,1),tau,MaxDim,Rtol,Atol,speed);

% Create a plot so we can see the results.
figure
plot(1:MaxDim, dim, 'o-k')
hold
plot(dE, dim(dE), 'or') % point where NN < 0.1 percent - overly restrictive?
xlabel('Embedding Dimension')
ylabel('% FNN')
grid on
   
% Wulf's algorithm
Fs = round(1/mean(diff(t))) % estimate sampling rate
period_samples = 48; % estimate period in samples (for evolve) (based on 4xTau) - there is also a power spectrum approach (see below)
tau = 12; % selected delay from AMI
m = 3; % selected embedding dimension from FNN
evolve = round(0.3*period_samples); % evolution between replacements (in samples). Fs is sampling rate (in Hz)

% Run Wulf's estimate for a small subset of the data to visualize the algorithm:
[out_w, lyew] = LyE_W_NK(y, Fs, tau, m, evolve,0);

% Plot segment LyE estimate and the final estimate
figure
plot(out_w(:,6))
xlabel('Segment')
ylabel('LyE estimate')
hold
line([0 max(out_w(:,1))],[lyew lyew],'color','r')

lyew

%% Another to get "typical" periodicity of the attractor 
% Using the power spectruum estimate
[Pxx,F] = periodogram(y(:,1),hamming(length(y)),length(y),60);

figure
plot(F,Pxx)
xlabel('Hz')
ylabel('Power (original units^2/Hz)')
title('Modified Periodogram Power Spectral Density Estimate')
