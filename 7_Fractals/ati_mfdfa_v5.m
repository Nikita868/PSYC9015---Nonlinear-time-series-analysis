function [generalizedHurst,multiSpectrum,MF_Width,MF_Dominant]  = ati_mfdfa_v5(z,qmin,qmax,npoints)

%% Inputs
% z - data (should be fGn type)
% qmin, qmax - minimal and maximal q-order exponents
% npoints - number of steps between qmin and qmax
% start_scale, end_scale - first and last scale for the linear fit to the
% fluctuation function

%% Outputs
% generalizedHurst - matrix contatining q values in the first row and
% q-generalized Hurst values in the second row

% multiSpectrum - matrix containing singularity strength D(h) [termed F(alpha) in
% Ihlen(2012)] in the first row. Holder exponent is in second row [termed
% alpha in Ihlen(2012)].

% MF_Width - simple estimate of multifractal spectrum width. Calculated as
% range of the Holder exponent values (max - min).

% MF_Dominant - estimate of average multifractal value (estimate of most
% dominant Holder exponent). Calculated at mean Holder exponent.

%% Notes
% Data should be integrated if the original looks more like an fGn signal
% and not integrated when it is more like fBm signal.

% Nikita Kuznetsov June 2018
% See Ihlen (2013) for an MF-DFA tutorial.

%% Example use
% data = ati_simulate_fGn(.6,2048,'fGn');
% [generalizedHurst,multiSpectrum,MF_Width,MF_Dominant]  = ati_mfdfa_v5(data,-3,3,10);


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
    
    
% Verify data is at least 1024 points in length
    if len_p2<10
        
        error('Series must be at least 1024 points in length.')
    
    else
    end
% % % Verify Z-score format, if not fix.
% %     if abs(mean(z))>=1e-6 || std(z,1)~=1
% %         
% %         mu=mean(z); sigma=std(z,1);
% %         z = (z - mu)./ sigma;
% %         disp('Data was normalized for analysis.')
% %         
% %     else
% %     end

data_orig= z;
    
%% Try to determine if it is a fGn or fBm
s_exp=spec_expon(z,(2^len_p2)*.25,.25);

if (.75 < s_exp) && (s_exp<=3)

  disp(['Signal might be fBm, spectral exponent= ' num2str(s_exp,2)]);
  pfile=input('Integrate time series? 0= no, 1=yes: ');
  if pfile==1
     z=cumsum(z); z=z-mean(z); % get the profile
  else end


elseif (-1 < s_exp) && (s_exp<=.75)

  disp(['Signal might be fGn, spectral exponent= ' num2str(s_exp,2)]);
   pfile=input('Integrate data? 0= no, 1=yes: ');
  if pfile==1
     z=cumsum(z); z=z-mean(z); % get the profile
  else end
else
  warning('Spectral analysis could not determine signal class (fBm or fGn)');
  disp('Assuming stationary fGn, results may be suspect!')

end
          
    
%% Specify the analysis, DFA or MFDFA and the different values of q
% type=input('Mf-DFA Analysis. Press return or 1 default q MFDFA [-10:10], or 2 to specify q: ');
  
% if isempty(type), type=1; end
%  
%     switch type
%         
%      case {1}
%             q=[-10, -8, -6, -5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 8, 10];
%       
%      case {2}
%             qmin=input('Specify minimum q (smallest moment): ');
%             qmax=input('Specify maximum q (largest moment): ');
%             npoints=input('Specify the how many q between qmin and qmax: ');
            q=linspace(qmin,qmax,npoints);
            q(q==0)=[];
%     end
 

figure
subplot(2,3,1:3)
plot(data_orig,'k.-')
ylabel('Amplitude, A(t)')
xlabel('Time, t')
xlim([0 length(data_orig)])
title('Original time series')
set(gca,'fontsize',12)

% this column vector stores the Hurst exponents for each q
 ho=zeros(1,size(q,2));
% this is the q loop, the analysis is redone at each value of q
for m=1:size(q,2)
             Fluctuation=zeros(len_p2-3,1); % this re-initalizes the mean variance vector
                                            % on each loop
             for j=2:(len_p2-2)
                % computes bin size for each scale
                bin_size=length(z)/pow2(j);
                k = fix(length(z)/bin_size);               
                segment_vars=var(detrend( reshape(z,bin_size,k) ),1);      
                %raises the variance to a power of q, computes mean, then
                %takes inverse of q
                Fluctuation(j-1,1)=mean(segment_vars.^(q(m)/2))^(1/q(m));
             end
            
             c=flipud(Fluctuation)';
             c(1:2)=[]; % eliminate 2 smallest scales
             
             p=polyfit(4:(len_p2-2),log2(c),1); % fit a line, slope is H.
             ho(1,m)=p(1);
             % plots as a function of scale, number of points IN a
             % bin, which is why the Fluctuation vecotor is Flipped
             % upside-down (i.e. matlabs flipud function)
             subplot(2,3,4)
             plot(2:(len_p2-2),log2(flipud(Fluctuation+eps)),'o-k'), hold on
             title('Generalized fluctuation function')
                          
             plot(4:(len_p2-2),polyval(p,4:(len_p2-2)),'r-')
             xlabel('log_2(bin size, n)')
             ylabel('log_2[F(n)^q]'), %axis([2 (len_p2-2) -5 10])
             set(gca,'fontsize',12)
end
h=ho;                   
hold off, drawnow


    % Compute various versions of the MF statistics  
    Tau=q.*h-1;
    %D=Tau./(q-1);
    dtau=diff(Tau)./diff(q);
    Tau1=dtau(1:length(dtau)-1);
    Tau2=dtau(2:length(dtau));
    Alpha=(Tau1+Tau2)./2;
    F=q(2:length(q)-1).*Alpha-Tau(2:length(q)-1);
    
    mf_width=range(Alpha);
    mf_m=median(Alpha);
    disp(['Range: ' num2str(mf_width,2) ' Median: ' num2str(mf_m,2)])
    % Make plots
    
    subplot(2,3,5)
    plot(q,h,'ok'), %axis([-10 10 .5 1.5])
    %hold on, plot(q(14),h(14),'rd')
    xlabel('Moment order, q');ylabel('H(q)')
    title('Generalized Hurst Exponent, H(q)')
    ylim([0 2])
    set(gca,'fontsize',12)
    
%     subplot(3,2,5)
%     plot(q, Tau,'ok'), %axis([-10 10 -15 15])
%     hold on
%     xlabel('q');ylabel('\tau(q)')
    
   % subplot(2,2,3)
  %  plot(q,D,'o'), %axis([-10 10 .5 1.5])
   % hold on
  %  xlabel('q');ylabel('Multifractal dimension D(q)')
    
    subplot(2,3,6)
     plot(Alpha,F,'o-k'), %axis([.5 1.75 .0 1.2])
     xlabel('Hölder exponent, h');ylabel('Singularity Strength, D(h)')
    title('Multifractal spectrum')
    ylim([0 1])
    xlim([0 2])
    set(gca,'fontsize',12)
    p2 = polyfit(Alpha,F,2)
    hold
    plot(Alpha,polyval(p2,Alpha),'r')
     
%      set(gcf,'Position',[11         484        1222         631])
     

%% Finish up, write outfile
   
Alpha=[Alpha,NaN,NaN];
F=[F,NaN,NaN];
% distMatrix=[q;h;Tau;Alpha;F];
multiSpectrum = [Alpha;F];
generalizedHurst = [q;h];


%% MF-spectrum width
MF_Width = range(multiSpectrum(1,:));
MF_Dominant = nanmean(multiSpectrum(1,:));


end

%Spectral routine to determine signal class
function [sexp]=spec_expon(z,nfreq,spec_portion)
%returns the slope of the powerspectrum
%function [pslp]=spec_slope(z,nfreq,spec_portion);
%default nfreq=128, default portion=1
if nargin == 1
	nfreq=128;
    spec_portion=1;  
else end

[P,F]=pwelch(z,triang(2*nfreq),nfreq,2*nfreq,1,'onesided');
%[P, F] = spectrum(z,2*nfreq,nfreq,triang(2*nfreq),1,.95,'none');
F(1)=[];
F(nfreq)=[];
P(1,:)=[];
%P(:,2)=[];
P(nfreq,:)=[];
P=log10(P);
F=log10(F);
[p,S]=polyfit(F(1:floor(length(F)*spec_portion)),P(1:floor(length(P)*spec_portion)),1);
sexp=-1*p(1,1);
end
     
     
