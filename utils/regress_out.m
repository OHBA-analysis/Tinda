function data_corrected = regress_out(data, regressors)
if ~exist('regressors', 'var') || isempty(regressors)
    warning('There were not regressors. data_out = data_in')
    data_corrected = data;
else
    if ~all(regressors(:,1)==1)
        regressors = [ones(size(regressors,1),1), regressors];
    end
    beta = regressors\data; % B = X\Y
    model = regressors * beta;
    data_corrected =  data - model;
end