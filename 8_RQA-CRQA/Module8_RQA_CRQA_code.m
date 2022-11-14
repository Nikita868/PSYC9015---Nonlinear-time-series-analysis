%% This code accompanies the powerpoint slides on Recurrence quantification

%% 1. Example 1: Application of RQA to thoracic chest expansion measurements during breathing

% Load data
load('breathing_spont_YA.mat')

% Sampling rate is 50 Hz and duration of recording is 30 seconds

% Plot the data - always plot the data
plot(t,data)
ylabel('Time (s)')
xlabel('Chest expansion (mm)')
legend('Thoracic','Abdominal')

% Find Time delay
L = 300;
[tau_s, ami_s] = AMI_Stergiou(data(:,1),L);

% Plot AMI function
figure
plot(ami_s(:,1),ami_s(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI';'Stergiou'})
grid minor

% Find embedding dimension
MaxDim = 12;
tau = 64; % this value is found from the AMI function above.

Rtol = 15;
Atol = 2;
speed = 0;
[dE,dim] = FNN(data(:,1),tau,MaxDim,Rtol,Atol,speed);

% FNN plot
figure
plot(1:MaxDim, dim, 'o-k')
hold
plot(dE, dim(dE), 'or') % point where NN < 0.1 percent - overly restrictive?
xlabel('Embedding Dimension')
ylabel('% FNN')
grid on

% Run RQA - first set RQA parameters
% This RQA code is a one-stop-shot with many parameter choices. These are
% the most generic you could choose. If you are using different data the
% SETVALUE parameter may need to be tweeked. Here I choose a value (5 deg)
% that produces an interesting RQA plot.
TYPE = 'RQA';
EMB = 3;
DEL = tau;
ZSCORE = 0;
NORM = 'EUC'
SETPARA = 'radius'
SETVALUE = .1
PLOTOPTION = 1;

% Run RQA on thoracic data
[RP,results] = RQA(data(:,1), TYPE, EMB, DEL, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);

% Compare to abdominal data

%% 2. Compare to faster breathing trial
clear

% Load data
load('breathing_paced_YA.mat')

% Sampling rate is 50 Hz and duration of recording is 30 seconds

% Plot the data - always plot the data
plot(t,data)
ylabel('Time (s)')
xlabel('Chest expansion (mm)')
legend('Thoracic','Abdominal')



%% 3. Categorical RQA

% import data
load('lyricsA.mat')

% Run RQA - first set RQA parameters
% This RQA code is a one-stop-shot with many parameter choices. These are
% the most generic you could choose. If you are using different data the
% SETVALUE parameter may need to be tweeked. Here I choose a value (5 deg)
% that produces an interesting RQA plot.
TYPE = 'RQA';
EMB = 1;
DEL = 1;
ZSCORE = 0;
NORM = 'NON'
SETPARA = 'radius'
SETVALUE = .005
PLOTOPTION = 1;
[RP,results] = RQA(lyricsA, TYPE, EMB, DEL, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);



%% 4. Cross-recurrence (needs update - 11/14/22)
clear
clc

% Here we are going to load joint angle data and run cRQA on two of the 
% time series. We will have to get the state space parameters from both and
% choose which to use. Some of the other parameters will also be changed
% from above.

load('.......')

% We'll only use 1000 data point to save processing time. you should step
% up what data lengths you run to get an idea of the processing time.
data1 = z(1:1000);
data2 = y(1:1000);
data3 = x(1:1000);
% This data was sampled at 60 Hz.
f = 60;

% AMI needs to be run to find a proper time lag.
L1 = f;
[tau1, tau_v1] = AMI_Stergiou(data1, L1);

L2 = f;
[tau2, tau_v2] = AMI_Stergiou(data2, L2);
% FNN needs to be run to find a proper embedding dimension.
MaxDim = 12;
Rtol = 15;
Atol = 2;
speed = 1;
[dE1, dim1] = FNN(data1, tau1(1), MaxDim, Rtol, Atol, speed);
[dE2, dim2] = FNN(data2, tau2(1), MaxDim, Rtol, Atol, speed);

% This data happens to get the same embedding dimension of 4. For the time
% lag will use an integer multiple of the two.
dim = dE1;
tau = round(tau1(1)*tau2(1));

% This is similar to above except the normalization and radius value were
% changed. It is important to concider what these mean when performing RQA.
% Here I choose a value (5 deg) that produces an interesting RQA plot.
TYPE = 'mdRQA';
EMB = dim;
DEL = tau;
ZSCORE = 0;
NORM = 'EUC';
SETPARA = 'radius';
SETVALUE = 0.5;
PLOTOPTION = 1;
results = RQA([data1,data2,data3], TYPE, EMB, DEL, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);





