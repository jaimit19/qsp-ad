function samples =morrisFdSampling(problem, N, delta)
arguments
    problem;
    N;
    delta =1e-3; % optional input parameter or vector 
end
% morrisFdSampling generates samples using lhs scheme to pefrom 
% Morris Finite Different based SA analysis.
% The samples are scaled by upper and lower bounds.
% Inputs
%   problem - struct
        % problem has to have the following fields:
        %    lowerBound - desired lower bound of each of the parameters
        %    upperBound - desired upper bound of each of the parameters
        %    num_vars - number of parameters 
        %    names - optionally name of the parameters 
%   N - int The total number of returned samples is N * num_vars
%   delta -  step size for finite difference, set to default value of 1 / 20
% Outputs
%   samples - The samples for morris finite different gsa analysis
% Last updated July 29 2021 by Jaimit Parikh


if ~all(isfield(problem, ["lowerBound", "upperBound", "num_vars"]))
    error("The necessary fields in the problem structure are not included."...
        +" The required field names are lowerBound, upperBound, and num_vars")
end

import sa.latinHypercubeSampling
iniSamples = latinHypercubeSampling(problem, N);

nPars = problem.num_vars; 

samplesN = repelem(iniSamples, nPars+1, 1);
deltaM =[zeros(1, nPars); delta.*diag(ones(1, nPars))];
deltaM = repmat(deltaM, N, 1);


samples = samplesN + deltaM ; 

end