function output = ati_pspc(z,window_size,spec_portion)
%Compute power spectrum using Welch's method See
%Holden J. G. (2005).  Gauging the fractal dimension of response times from cognitive tasks.
%    In M. A. Riley & G.  C. Van Orden (Eds.), Contemporary nonlinear methods for behavioral
%    scientists: A webbook tutorial, 267-318. At http://www.nsf.gov/sbe/bcs/pac/nmbs/nmbs.pdf
%Jay Holden, 7/2006

% Inputs
% 1) window_size. This referst to the window size for FFT - in segments as an integer power of 2.
% Nfreq=(segment length/2+1),Max is N/2,Nyquist and DC are deleted for n-1
% Default is 1/8*(total length), shift is Nfreq. Asuming N=1024, averages of 7 segments are output.

% 2) spec_portion - What percentage of low frequencies to fit? (e.g. 50%)

disp('Specify Number of frequencies to estimate')
disp('in segments as an integer power of 2.')
disp('Nfreq=(segment length/2+1),Max is N/2,')
disp('Nyquist and DC are deleted for n-1')
disp('Default is 1/8*(total length), shift is Nfreq.')
disp('Asuming N=1024, averages of 7 segments are output.')

nfreq=window_size;
spec_portion=spec_portion/100


data = z;

%% Error Checking, Verify single column of data
    [zr,zc] = size(z);
    if zr > 1 && zc > 1
         error('Matrix Input: Assumes a single column of data.')
         
    elseif zr==1    

        z = z';
        warning('Row Input: Assumes a single column of data.')
        
    else   
    end

% Verify data is an integer power of 2 in length
    len_p2=log2(length(z));  % find out how long it is
    
    if mod(len_p2,1)~=0 %Integer test, make sure no decimal remainder after dividing by 1

        error('Series must be an integer power of 2 in length.')
        
    else
    end
% Verify data is at least 64 points in length
    if len_p2<6
        
        error('Series must be at least 64 points in length.')
    
    else
    end
% Verify Z-score format, if not fix.
    if abs(mean(z))>=1e-6 || std(z,1)~=1
        
        mu=mean(z); sigma=std(z,1);
        z = (z - mu)./ sigma;
        disp('Data was normalized for analysis.')
        
    else
    end
    
%get the number of frequencies

    if nfreq==[]
        
        nfreq=.125*length(z);
        
    elseif mod(log2(nfreq),1)~=0
        
        error('Nfreq must be an integer power of 2.')
        
    else
    end
    
%% Begin Power Spectrum
%[Pxx,w] = pwelch(x,window,noverlap,nfft)

[P,F]=pwelch(z,triang(2*nfreq),nfreq,2*nfreq,1,'onesided');
F(1)=[];
F(nfreq)=[];
P(1,:)=[];
P(nfreq,:)=[];

%% Plot Figure
figure
subplot(121)
plot(F,P,'.-')
ylabel('Power')
xlabel('Frequency')
subplot(122)
plot(log10(F),log10(P),'.-')
ylabel('Log_1_0(Power)')
xlabel('Log_1_0(Frequency)')
% p = polyfit(log10(F),log10(P),1);
[p]=polyfit(log10(F(1:floor(length(F)*spec_portion))),log10(P(1:floor(length(P)*spec_portion))),1);

hold
plot(log10(F(1:floor(length(F)*spec_portion))),polyval(p,log10(F(1:floor(length(F)*spec_portion)))),'r','linewidth',2)

beta_r = round(p(1)*100)/100;
text(min(log10(F)),mean(log10(P)),['\beta = ' num2str(beta_r)])

if -p(1) > 1;
title(['Series fBm. H = ' round(num2str((-p(1)-1)/2)*100)/100])
else
title(['Series fGn. H = ' round(num2str((-p(1)+1)/2)*100)/100])    
end

% set(gcf,'Position',[12   737   992   379])

%% Finish up
 %put in the outfile
 P=log10(P);
 F=log10(F);
 output=[F P]';

 
 
 