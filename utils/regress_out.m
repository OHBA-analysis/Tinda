function data_corrected = regress_out(data, regressors)

beta = regressors\data; % B = X\Y
model = regressors * beta;
data_corrected =  data - model;