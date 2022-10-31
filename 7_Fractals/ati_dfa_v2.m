function [output,Alpha] = ati_dfa_v2(data,minscale,maxscale)
%% Function for performing detrended fluctuation analysis
% Copyrighted by Jing Hu and Jianbo Gao, ECE department, University of Florida
% this function was re-adapted for ATI 2015 by Nikita Kuznetsov. May 2015.


%% Inputs
% data - either single row or single column
% minscale - smallest scale for fitting (in log 2 coordinates) - e.g. 2
% means 4 samples
% maxscale - largest scale for fitting.
% integrate - perform integration? 1 - yes; 0 - no


%% Outputs
% output - Fluctuation function vs. scale plot
% slope - DFA Alpha, equivalent to slope

%% Example use
% data = ati_simulate_fGn(.4,2048,'fGn');
% [output,alpha] = ati_dfa_v2(data,4,9);

% Copyrighted by Jing Hu and Jianbo Gao, ECE department, University of Florida
% This function was re-adapted for Nonlinear Methods ATI 2015 by Nikita Kuznetsov. May 2015.

%% make sure data are in column format.
if size(data, 1) > size(data, 2)
    data = data';
end

data_or = data; % keep original data for plotting
data = data - mean(data);

Len = length(data);

% minscale = 2;
% maxscale = floor(log2(Len));
step = 0.5;

%% Integrate data (only do if fGn-like) - we integrate in the data input only if needed
Y = data;
    
% divide the integrated time series into boxes of equal length
k = 1;
for i = minscale : step : maxscale
    x = Y;
    l = floor(2^i);
    
    b_block = floor(Len/l);
    clear z;
    
    for j = 1 : b_block
        zx = 1 : l;
        zy = x((j - 1)*l + 1 : j*l);
        sx = sum(zx);
        sy = sum(zy);
        sxy = sum(zx.*zy);
        sx2 = sum(zx.^2);
        x_bar = sx/l;
        y_bar = sy/l;
        slope = (l*sxy-sx*sy)/(l*sx2-sx*sx);
        zz = slope*(zx - x_bar) + y_bar;
        z((j - 1)*l + 1 : j*l) = (x((j - 1)*l + 1 : j*l) - zz).^2;
    end
    dfa_data(k) = (sum(z)/length(z))^(1/2);
    k = k + 1;
end

result(1, :) = minscale : step : maxscale;
result(2, :) = log2(dfa_data);

%% Plot Figure
figure
subplot(121)
plot(data_or,'.-k')
ylabel('Amplitude, A(t)')
xlabel('Time, t')
xlim([0 length(data)])

subplot(122)
plot(result(1,:),result(2,:),'ko-','linewidth',2)
ylabel('Log_2(F(n))')
xlabel('Log_2(Window size, n)')
% set(gcf,'Position',[12   737   992   379])
axis square
grid on

p = polyfit(result(1,:),result(2,:),1);
hold
plot(result(1,:),polyval(p,result(1,:)),'r','linewidth',2)

alpha_r = round(p(1)*100)/100;
if p(1) > 1
    H = round((p(1)-1)*100)/100;
    info = 'Series fBm.';
elseif p(1) < 1
    H = alpha_r;
    info = 'Series fGn.';
end
title(['\alpha = ' num2str(alpha_r) ': ' info ' H = ' num2str(H)])

output=[result(1,:) result(2,:)];
Alpha = p(1);


