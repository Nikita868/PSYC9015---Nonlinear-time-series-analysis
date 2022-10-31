% Fractal module sample code:

%% Simulate fGn and fBm signals

H = 0.5; % Hurst exponent
sim_data = ati_simulate_fGn(H,1024,'fGn');

% plot simulated data
figure
plot(sim_data,'k.-','linewidth',1)
xlabel('Time, t')
ylabel('Amplitude, Y(t)')
xlim([0 length(f)])
% ylim([-5 5])
box off
% set(gcf,'Position',[11   758   960   356])
title(['Simulated fGn. H = ' num2str(round(H*100)/100)])
set(gca,'PlotBoxAspectRatio',[5,1,1])

%% Spectral analysis on simulated signals:
%Compute power spectrum using Welch's method See
%Holden J. G. (2005).  Gauging the fractal dimension of response times from cognitive tasks.

% 1) window_size. This referst to the window size for FFT - in segments as an integer power of 2.
% Nfreq=(segment length/2+1),Max is N/2,Nyquist and DC are deleted for n-1
% Default is 1/8*(total length), shift is Nfreq. Asuming N=1024, averages of 7 segments are output.

% 2) spec_portion - What percentage of low frequencies to fit? (e.g. 50%)
[output] = ati_pspc(sim_data,256,50);


%% SDA application to simulated data

% Note that the input data to SDA needs have length as a power of 2...
output = ati_sda(sim_data);

% Plot results
figure
subplot(121)
plot(sim_data,'k.-')
ylabel('Amplitude, Y(t)')
xlabel('Time, t')
xlim([0 length(sim_data)])
 
subplot(122)
plot(output(:,1),output(:,2),'ko-','linewidth',2)
hold
xlabel('Log_2(Bin Size)')
ylabel('Log_2(SD Dispersion)')
p = polyfit(output(1:end-2,1),output(1:end-2,2),1);
plot(output(1:end-2,1),polyval(p,output(1:end-2,1)),'r','linewidth',2)
set(gcf,'Position',[1         650        1122         466])

p(1)=round(p(1)*100)/100;
title(['SDA Slope = ' num2str(p(1)) ' H = ' num2str(p(1)+1) ' [H = Slope + 1]'])

%% Detrended Fluctuation Analysis

% Let's test the DFA algorithm on our simulated signal: 
% input parameters:
% data 
% minscale - smallest scale for fitting (in log 2 coordinates) - e.g. 2 means 4 samples
% maxscale - largest scale for fitting.
[output,Alpha] = ati_dfa_v2(integrate(sim_data),2,9);

% Note that you only need to use the integrate function on your data if the
% signal is like an fGn. If it's like fBm, then no integration is necessary

%% Adaptive Fractal Analysis
% Visualize the continuous trend fit using AFA
[detrended_data, record_y,record_x] = detrending_method_nk(integrate(sim_data), 2^6+1, 1);

% AFA application
[result1,p1_slope, intercept] = AFA_function(integrate(sim_data),2,9)

% Note that you only need to use the integrate function on your data if the
% signal is like an fGn. If it's like fBm, then no integration is necessary

%% Multifractal DFA
% Simulate a slightly longer fGn process with H = 0.6
H = 0.6;
sim_data = ati_simulate_fGn(.6,2048,'fGn');

% Run MF-DFA:

%% Inputs
% z - data (should be fGn type)
% qmin, qmax - minimal and maximal q-order exponents
% npoints - number of steps between qmin and qmax
% start_scale, end_scale - first and last scale for the linear fit to the
% fluctuation function
[generalizedHurst,multiSpectrum,MF_Width,MF_Dominant]  = ati_mfdfa_v5(sim_data,-3,3,10);



