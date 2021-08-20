function S = morrisFdAnalysis(X, Y)
arguments
    X; %  n x k ... n input samples with k parameters
    Y; % n x m ... n samples with m output variable of interest  
end
% Morris finited difference based analysis to estimate the sensitivity of
% response variable Y to parameters X
% Inputs
    % X -   n x k ... n input samples with k parameters
    % Y - n x m ... n samples with m output variable of interest  
% Outputs
    % S - strucutre with sensitivity indices
    %   S.mustar - 
    %   S.sigma -


% Last updated July 29th 2021 by Jaimit Parikh
k = size(X, 2);  %number of parameters
r = size(X, 1) / (k + 1);
yStar = Y(1:k+1:end, :);
mu = zeros( k, size(Y,2));
mus = zeros(k, size(Y,2));
sigma= zeros(k, size(Y,2));
variation = max(X) - min(X);
for ii = 1:k
    delta = X(ii+1, ii) - X(1, ii);
    disp(delta);
    yA = Y(ii+1:k+1:end, :);
    fd = variation(ii) * (yA - yStar) / delta ; 
    mu(ii, :) = sum(fd) / r;
    mus(ii, :) = sum(abs(fd)) / r ;
    sigma(ii, :)  = sum((fd - mu(ii)).^2) / (r - 1); 
end

  S.mustar = mus;
  S.sigma = sigma;

end