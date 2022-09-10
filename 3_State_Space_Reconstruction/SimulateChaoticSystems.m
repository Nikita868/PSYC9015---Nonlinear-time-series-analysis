%% Simulate nonlinear systems
% Wrapper script to simulate several nonlinear dynamical systems with
% chaotic dynamics

% Nikita Kuznetsov 9/7/2022 - kuznetna@ucmail.uc.edu

%% 1.Simulate nonlinear systems
% Below is the code to simulate several nonlinear dynamical systems with
% chaotic dynamics. It uses the the ChaosLibrary MATLAB m file uses ODE
% solvers to calculate various chaotic attractors. You can change all
% components of the chaotic attractor, including the time scale of
% simulation, the initial conditions and the internal parameters of the
% attractor equations.

% inputs  - S, string name of chaotic attractor
%         - t, time, either of the form [to,tf] or [to:f:tf]
%         - IC, initial conditions
%         - p, coefficients of the system of differential equations used.
% outputs - t, time
%         - y, column oriented time series

% Simulate the Lorenz system (see the ChaosLibary.m for the explanation of
% the inputs and other possible systems to simulate.
[t,y] = ChaosLibrary('Lorenz',[0,100],[0.1 -0.01 9],[10 28 8/3]);

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

% Plot the Lorenz system state space
figure
plot3(y(:,1),y(:,2),y(:,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
axis square
hold
view(55,31)

% Animated plot of the trajectory
figure
grid on
comet3_2(y(:,1),y(:,2),y(:,3))
% comet3_2.m is equivalent to standard matalab comet3, but implements a
% pause in the drawing such that the trajectory evolves slower.

%**************************************************************************

% Simulate Rossler system
[t,y] = ChaosLibrary('Rossler',[0,100],[-9 0 0],[0.2 0.2 5.7]);

% plot individual variables of Rossler
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

% plot state space of Rossler
figure
plot3(y(:,1),y(:,2),y(:,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
axis square
hold

% Animated plot of the trajectory
figure
grid on
comet3_2(y(:,1),y(:,2),y(:,3))

%% 2. Average Mutual Information and Choosing time delay for embedding

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

% This will run the AMI_Thomas method (kernel density) with a max lag of
% 150.
[tau_t, ami_t] = AMI_Thomas(y(:,1), L);

% Plot AMI function 
figure
plot(ami_t(:,1),ami_t(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI';'Thomas'})
grid minor

% Compare histogram and kernel density function methods
figure
plot(ami_s(:,1),ami_s(:,2), 'b')
hold
plot(ami_t(:,1),ami_t(:,2), 'k')
xlabel('Time Lag')
ylabel({'AMI'})
legend('Stergiou (histogram)','Thomas (density)')
grid on

%% Automated peak detection of the first minimum in the AMI function
% Estimate location of the first minimum - this would not always work, but
% it's a start. Further inputs could be added to the findpeaks function to
% finetune performance for specific datasets.
[PKS,LOCS] = findpeaks(-ami_t(:,2),'NPeaks',1);

figure
plot(ami_t(:,1),ami_t(:,2), 'k')
hold
plot(LOCS,ami_t(LOCS,2),'*r')
xlabel('Time Lag')
ylabel({'AMI';'Thomas'})
grid minor

%% 3. Finding the embedding dimension using the FNN algorithm
% The maximum dimension is largely determined through experience. The
% method is pretty quick so a higher 12 doesn't add much processing time.
% The other parameters come from the publication. 
% - Reference:   "Determining embedding dimension for phase-space
%                 reconstruction using a geometrical construction",
%                 M. B. Kennel, R. Brown, and H.D.I. Abarbanel,
%                 Physical Review A, Vol 45, No 6, 15 March 1992,
%                 pp 3403-3411.
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


%% Reconstruct the attractor using embedding and tau
y_rec = embed(y(:,1),3,12)';

% plot state space of Rossler
figure
subplot(121)
plot3(y(:,1),y(:,2),y(:,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
axis square
title('Rossler')

subplot(122)
plot3(y_rec(:,1),y_rec(:,2),y_rec(:,3))
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
axis square
title('Reconstructed Rossler')




