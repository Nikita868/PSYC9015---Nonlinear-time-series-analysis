function data = rm_outliers(data,sd_value)
% Remove outliers above and below mean +-SD 
% Inputs:
% data - data in column format
% sd_value - number of SDs for outlier detection

mu = mean(data);
sigma = std(data);
data(abs(data - mu) > sd_value*sigma)=[];

end