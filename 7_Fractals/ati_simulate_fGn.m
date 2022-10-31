function f = ati_simulate_fGn(H,N,signalClass)
% This function simulates fGn time series using Davies and Harte's (1997) algorithm. Steps
% in the algorighm were presented in Caccia et al. (1996). The
% implementation of steps 3 and 4 is borrowed from the ffgn function by 
% Yingchun Zhou (Jasmine), zhouyc@math.bu.edu and Stilian Stoev, sstoev@umich.edu)

%% Inputs:
% H - Hurst exponent value. Should be between 0 and 1. 
% N - Desired length of time series. 
% signalClass - Desired class of the signal ('fBm' or 'fGn')
%% Outputs:
% f - Time series of lenght N.


%% File Saving prompts here:
% [OutFileName,OutPathName] = uiputfile('*.txt','Specify an output file name');
% 
% if isequal(OutFileName,0) || isequal(OutPathName,0)
%    %quit if no outfile is specified
%     error('Canceled')
% end
    
% %% User input
% H = input('Specify H coefficient (between 0 and 1): ');
% if abs(H) >= 1
%     error('H should be between 0 and 1')
% end

% N = input('Desired time series length (has to be power of 2): ');
% if mod(log2(N),1)~=0 %Integer test, make sure no decimal remainder after dividing by 1
%    disp('Series must be an integer power of 2 in length.')
%    N=2^nextpow2(N);
%    disp(['N set to ' num2str(N) ])
% end

%% May 2015 Nikita Kuznetsov


%% Specify the desired autocorrelation
M = 2*N;
n = 1;
sigma = 1;
tau = 0:M/2;
gamma = (sigma^2/2).*(abs(tau+1).^(2*H)-2.*abs(tau).^(2*H)+abs(tau-1).^(2*H));
g = [gamma(1:end) fliplr(gamma(2:end-1))];

%% Step 1 (Calculate exact spectral power expected for the theoretical autocorrelation function)
S = real(fft(g));

%% Step 2 (Check if all spectral coefficients are greater than 0)
if min(S)<0
error(' Some of the Sj are negative!');
end
S = abs(S);

%% Step 3 (Calculate randomized spectral amplitudes) 
z(:,1)=sqrt(2)*randn(n,1);
y(:,1)=z(:,1);

z(:,N+1)=sqrt(2)*randn(n,1);
y(:,N+1)=z(:,N+1);

a=randn(n,N-1);
b=randn(n,N-1);

z1=a+b*1i;
z(:,2:N)=z1;
y1=z1;

y(:,2:N)=y1;
y(:,N+2:2*N)=conj(y(:,N:-1:2));

y = y.*(ones(n,1)*sqrt(S));

%% Step 4 (use
f = real(fft(y')')/sqrt(4*N);
f = f(:,1:N);
% plot(f)

% Random number generator for H = 0.5
if H == .5
    f = randn(N,1);
end


%% Convert simulated fGn into fBm if needed
if  strcmp(signalClass,'fBm') == 1  
    % Cumulative sum of the data. This performs the integration.
    f = cumsum(f);
else
    f = f;
end

   
end

