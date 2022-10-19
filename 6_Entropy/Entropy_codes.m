%% Entropy chapter
% Here are some examples of code to calculate variour types of entropies:
%   Sample Entropy
%   Permutation Entropy
%   Symbolic Entropy
%   Multiscale entropy
%   Refined Composite Multiscale Entropy.

%% Sample Entropy
% Source code for SampEn is from PhysioNet: https://physionet.org/content/sampen/1.0.0/

% A test signal is included in the text file c/sampentest.txt
% and the Matlab file matlab/*/sampentest.mat. The sample entropy
% calculations are implemented both in Matlab and a command-line
% executable obtained from C source code. The following Matlab
% session illustrates how to use each method and they give essentially
% the same result. Note that the first line of output is Sampen(0,r,N)
% which corresponds to m=0 and can be interpreted as the negative
% logarithm of the probability of a match of length 1.

% load test signal
load sampentest.mat

% run SampEn: 
%Input Parameters:
%sflag    flag to standardize signal(default yes/sflag=1) 
%cflag    flag to use fast C code (default yes/cflag=1) 
%vflag    flag to calculate standard errors (default no/vflag=0) 
%
%Output Parameters:
%se standard error estimates for m=0,1,...,M-1
%A number of matches for m=1,...,M
%B number of matches for m=0,...,M-1
%  (excluding last point in Matlab version)
%
% Open the sampen() function to see the descriptions for the
% input and output parameters
[e,se,A,B]=sampen(z,5,.2,1,0,1);

% overall outoput - it should be identical to the results in https://physionet.org/content/sampen/1.0.0/00README
e

% Sample entropy estimate for template length m=2
e(3)

% Sample entropy estimate standard error for template length m=2
se(3)

% Number of matches of m=2 and m=3
A(3) % m+1 = 3 matches
B(3) % m = 2 matches

%% Another SampEn calculation 
% This one is based on artificial data in 
% Kuznetsov, N., Bonnette, S., & Riley, M. A. (2013). Nonlinear time series
% methods for analyzing behavioural sequences. In Complex systems in sport
% (pp. 111-130). Routledge.

load SampleData.mat

figure
subplot(211)
plot(x(:,1),'.-')
subplot(212)
plot(x(:,2),'.-')

% Analyze the first time series (column 1 of x)
[e,se,A,B]=sampen(x(:,1),5,1,0,0,1); % note that we do not standardize the time series here (sflag = 0)

% Sample entropy estimate for template length m=2 and its standard error
e(3) % see the chapter for the expected value
se(3) % note that the standard error is very high! the data lenght is too short

% Analyze the second time series (column 1 of x)
[e,se,A,B]=sampen(x(:,2),5,1,0,0,1); % note that we do not standardize the time series here (sflag = 0)

% Sample entropy estimate for template length m=2 and its standard error
e(3) % see the chapter for the expected value
se(3) % note that the standard error is very high! the data lenght is too short


%% SampEn on stride data
load step_time.mat

figure
plot(x)
ylabel('Time (s)')
xlabel('Step (N)')

% Sampen with m = 2, r = 0.2 of standard deviation of the data
[e,se,A,B]=sampen(x,5,.2,1,0,1); 

% Sample entropy estimate for template length m=2 and its standard error
e(3) 
se(3) 

%% Permutation entropy

% These data are from the slides
data = [5 8 2 1 9 7 3 1 8];

% inputs -  data: 1-D array of data being analyzed
%           m: embedding dimension (order of permutation entropy) 
%           tau: time delay
% outputs - permuEnt: value calculated using a log base of 2
%           hist: number of occurences for each permutation order
m = 2;
tau = 1;
[permEnt, hist] = Ent_Permu(data, m, tau);

%% Symbolic entropy

load SymbolicEnt_data.mat

figure
plot(data,'o-')

% Convert the time series into a binary symbol series based on a threshold
% value
data_bin = data>0.5

% [ SymEnt ] = Ent_Symbolic20180320( X, L )
% symbolicEnt Calculates the Symbolic Entropy for given data.
% Input -   X: 1-Dimensional binary array of data
%           L: Word length
% Output -  NCSE: Normalized Corrected Shannon Entropy

NCSE = Ent_Symbolic(data_bin,3)


%% Multiscale entropy

% simulate white noise
x = randn(1000, 1);

% - This code finds the Refined Composite Multiscale Sample Entropy,
%   Composite Multiscale Entropy, Multiscale Entropy, Multiscale Fuzzy
%   Entropy and Generalized Multiscale Entropy of a data series using the
%   methods described by - Wu, Shuen-De, et al. 2014. "Analysis of complex
%   time series using refined composite multiscale entropy." Physics
%   Letters A. 378, 1369-1374.
% - Each of these methods calculates entropy at different scales. These
%   scales range from 1 to tau in increments of 1.
% - The Complexity Index (CI) is not calculated by this code. Because the scales
%   are incremented by 1 the C is the summation of all the elements in each
%   array. For example the CI of MSE would be sum(MSE).

max_scale = 10;
m = 2;
[ RCMSE, CMSE, MSE, MSFE, GMSE ] = Ent_MS_Plus(x,max_scale,m,.2);

% Plot MSE result
plot(1:max_scale,RCMSE)
xlabel('Scale')
ylabel('Sample Entropy')



