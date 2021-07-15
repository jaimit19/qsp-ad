function S = regAnalysis(X, Y, regType, linearArgs)
arguments
    X; %  n x k ... n input samples with k parameters
    Y; % n x 1 ... nsamples with 1 output variable of interest ...  Currently 
    regType = "linear"; % available options are
    linearArgs.Intercept = true;
end

% Linear regression based analysis to estimate the sensitivity of
% the parameters X to the response variable Y 
% Inputs:
%   X - nxk matrix .. n input samples with k parameters
%   Y - n x1 matrix .. n samples with 1 output variable of interest
%   options - string for the kind of  correlation analysis to perform. default =
%   "partialcorr", other option include "corrcoef";
% Ouputs:
 %  S - output strucutre of sensitivity analysis
 %      S.S1 - normalized values of the sensivitiy indices
% Code by Jaimit Parikh

regType = string(regType);

X = normalize(X);
Y = normalize(Y);


switch regType
    case "linear"
        linearArgsCell = namedargs2cell(linearArgs);
        mdl = fitlm(X, Y, linearArgsCell{:});
        coefficients = mdl.Coefficients;
        coefs  = coefficients.Estimate;
        %coefs = (coefficients.Estimate).^2;
        S.S1 = coefs(2:end);
        S.Rsquared = mdl.Rsquared;
        fprintf("R2 score is: %d \n",  mdl.Rsquared.Adjusted);
end

end