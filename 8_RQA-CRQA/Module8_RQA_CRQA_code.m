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

TYPE = 'RQA'; % a string indicating which type of RQA to run (i.e.
%                 'RQA', 'cRQA', 'jRQA', 'mdRQA'). The default value is
%                 TYPE = 'RQA'.
EMB = 3; % the number of embedding dimensions (i.e., EMB = 1 would
%                be no embedding via time-delayed surrogates, just using
%                the provided number of colums as dimensions. The default
%                value is EMB = 1.
DEL = tau; % the delay parameter used for time-delayed embedding (if
%                EMB > 1). The default value is DEL = 1.
ZSCORE = 0; %indicates, whether the data (i.e., the different
%                columns of DATA, being the different signals or
%                dimensions of a signal) should be z-scored before
%                performing MdRQA:
%                0 - no z-scoring of DATA
%                1 - z-score columns of DATA
%                The default value is ZSCORE = 0.
NORM = 'EUC'; %NORM, the type of norm by with the phase-space is normalized.
%                The following norms are available:
%                'EUC' - Mean Euclidean distance norm
%                'MAX' - Maximum distance norm
%                'MIN' - Minimum distance norm
%                'NON' - no normalization of phase-space
SETPARA = 'radius'; % the parameter which you would like to set a target
%                value for the recurrence plot (i.e. 'radius' or
%                'recurrence'). The default value is SETPARA = 'radius'.
SETVALUE = .3; %  sets the value of the selected parameter. If
%                SETVALUE = 1, then the radius will be set to 1 if SETPARA
%                = 'radius' or the radius will be adjusted until the
%                recurrence is equal to 1 if SETPARA = 'recurrence'. The
%                default value if SETPARA = 'radius' is 1. The default
%                value if SETPARA = 'recurrence' is 2.5.
PLOTOPTION = 1; %  is a boolean where if true the recurrence plot will
%                be created and displayed.

% Run RQA on thoracic data
[RP,results] = RQA(data(:,1), TYPE, DEL, EMB, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);

% Optimize radius
rad = 0.05:.05:.8;
clear REC_range

for i = 1:length(rad)
    [RP,results] = RQA(data(:,1), TYPE, DEL, EMB, ZSCORE, NORM, SETPARA, rad(i), 0);
    REC_range(i,1) = [results.REC];
    disp([num2str(i) ' out of ' num2str(length(rad))])
end

figure
loglog(rad,REC_range,'o-')
xlabel('log(Radius)')
ylabel('log(%REC)')
ax = gca;
ax.XMinorTick = 'on';
ax.XAxis.TickValues = rad;
xtickangle(90)
ax.YMinorTick = 'on';
ax.YAxis.TickValues = REC_range;
grid on
axis square

%% (now you try!) Run the RQA analysis on the abdominal data

clear

% Load data
load('breathing_paced_YA.mat')

% Sampling rate is 50 Hz and duration of recording is 30 seconds

% Plot the data - always plot the data
plot(t,data)
ylabel('Time (s)')
xlabel('Chest expansion (mm)')
legend('Thoracic','Abdominal')


% RQA code here
%
%
%
%


%% Compare to faster breathing trial
clear

% Load data
load('breathing_paced_YA.mat')

% Sampling rate is 50 Hz and duration of recording is 30 seconds

% Plot the data - always plot the data
plot(t,data)
ylabel('Time (s)')
xlabel('Chest expansion (mm)')
legend('Thoracic','Abdominal')

% RQA code here
%
%
%
%

%% 2. Categorical RQA

% import data
load('lyricsA.mat')


% take a look at the data - try to understand if the analysis at the level
% of words or letters/ plot the data


% Run RQA - first set RQA parameters
TYPE = 'RQA';
EMB = 1; % embedding 1 - no embedding
DEL = 1; % time delay 1 - no embedding
ZSCORE = 0;
NORM = 'NON' % do not normalize the attractor (time series here)
SETPARA = 'radius'
SETVALUE = .005 % some small value should do
PLOTOPTION = 1;

[RP,results] = RQA(lyricsA, TYPE, DEL, EMB, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);

%% 4. Cross-recurrence
clear
clc

% Here we are going to load joint angle data and run cRQA on two of the 
% time series. We will have to get the state space parameters from both and
% choose which to use. Some of the other parameters will also be changed
% from above.

load('coord0.mat')

% This data was sampled at 120 Hz.
f = 120;

% Plot the data
figure
plot(data(:,1))
hold
plot(data(:,2))
ylabel('Rod tip position (mm)')
xlabel('Time (s)')
legend('Left hand','Right hand')
set(gcf,'Position',[-627,1284.20000000000,1060,367.200000000000])


%% AMI and FNN - do this for each time series separately and note the values
% AMI needs to be run to find a proper time lag.
L = 300
[tau_s, ami_s] = AMI_Stergiou(data(:,1), L);
[tau_s2, ami_s2] = AMI_Stergiou(data(:,2), L);

% Plot AMI function
figure
plot(ami_s(:,1),ami_s(:,2), 'k')
hold
plot(ami_s2(:,1),ami_s2(:,2), 'r')

xlabel('Time Lag')
ylabel({'AMI';'Stergiou'})
legend('data1','data2')
grid minor

% Find embedding dimension
MaxDim = 12;
tau = 30; % this value is found from the AMI function above.

Rtol = 15;
Atol = 2;
speed = 0;
[dE,dim] = FNN(data(:,1),tau,MaxDim,Rtol,Atol,speed);
[dE2,dim2] = FNN(data(:,2),tau,MaxDim,Rtol,Atol,speed);

% FNN plot
figure
plot(1:MaxDim, dim, 'o-k')
hold
plot(dE, dim(dE), 'or') % point where NN < 0.1 percent - overly restrictive?
plot(1:MaxDim, dim2, 'o-r')
xlabel('Embedding Dimension')
ylabel('% FNN')
legend('data1','data2')
grid on

%RQA parameters
TYPE = 'cRQA';
EMB = 3;
DEL = tau;
ZSCORE = 0;
NORM = 'EUC';
SETPARA = 'radius';
SETVALUE = 0.2;
PLOTOPTION = 1;

% Run RQA
[RP,results] = RQA([data(:,1),data(:,2)], TYPE, EMB, DEL, ZSCORE, NORM, SETPARA, SETVALUE, PLOTOPTION);






