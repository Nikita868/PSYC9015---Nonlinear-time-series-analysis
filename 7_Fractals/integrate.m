function data_int = integrate(data)
data = data - mean(data);
data_int = zeros(size(data));
data_int(1) = data(1);
for i = 2 : length(data)
    data_int(i) = data_int(i - 1) + data(i);
end