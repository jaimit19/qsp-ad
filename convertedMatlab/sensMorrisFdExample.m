% sensMorrisFdExample - Example file by Jaimit Parikh to run the  Finite difference 
% based morris sensitivity analysis on the QSP AD model based on code
% shared by Ralph Smith


import qspADmodel.*; % Import the functions from the qspADmodel package
import sa.*; % Import the functions from the sa (sensitivity analysis) package


% 1. Get parameters of the AD model
parameters = adModelPars();

% 2. Sample using Latin hypercube sampling the parameter values and running
% simulations for the samples

% Sampling for the desired parameters
parsName = ["k1",  "k2", "b1", "d1", "k5", "k6"];
nParameters = length(parsName);
problem.lowerBound = [2.6,  0.24,  6.3e-5,  0.15, 10, 25];
problem.upperBound = [10.4, 0.96, 2.52e-4, 0.6,  40, 100];
problem.num_vars = nParameters;
problem.names = parsName;
n = 300; % number of samples
delta = 1e-7; %(problem.upperBound - problem.lowerBound) / 1000;
samples = morrisFdSampling(problem,  n, delta);
N = size(samples, 1); 

% Simulating AD model for the desired parameter set and extracting SB and
% Th1 state variables
de = adDrugEffects();
IC = getInitialConditions();
t0 = 0; tfinal = 50;
SB = zeros(1, N); Th1 = zeros(1, N);
pN = cell(1,N);
for ii = 1:N
    pN{ii} = updateADPars(parameters, parsName, samples(ii, :));
    [T,Y] = ode15s(@(t, y)odefun(t, y, pN{ii}, de),...
        [t0 tfinal], IC);
    SB(ii) = Y(end, 1);
    Th1(ii) = Y(end, 3);
end

% Scatter Plot of SB and Th1 vs selected parameters
plotScatter(samples, SB, Th1, parsName, './figures/scatter.png');

Y = [SB', Th1'];
S = morrisFdAnalysis(samples, Y);

S1 = S.mustar ./ sum(S.mustar);

plotSAresults(S1(:, 1), S1(:, 2), parsName, 'temp.png', 'Morris FD');


function plotScatter(samples, SB, Th1,  parsName, fname)
f = figure('DefaultAxesFontSize', 14, 'Position', [40, 40, 1400, 800]);
for ii = 1:6
    subplot(2, 6, ii);
    plot(samples(:, ii), SB, 'ko'); xlabel(parsName(ii)); ylabel('SB');
    
    subplot(2, 6, ii+6);
    plot(samples(:, ii), Th1, 'ko'); xlabel(parsName(ii)); ylabel('Th1');
end
exportgraphics(f, fname, 'resolution', 300);
end

function p = updateADPars(p, parsName, parsValue)

for ii = 1:length(parsName)
    p.(parsName(ii)) = parsValue(ii);
end

end

function plotSAresults(SB, Th1, parsName, fname, tname)
f = figure('DefaultAxesFontSize', 14, 'Position', [40, 40, 1400, 800]);
subplot(1, 3, 1);
bar(SB); xticklabels(parsName); ylabel('PRCC');
title('SB')
subplot(1, 3, 2);
bar(Th1); xticklabels(parsName); ylabel('PRCC');
title('Th1');
subplot(1, 3, 3);
heatmap(["SB", "Th1"], parsName, [SB, Th1], 'Colormap', hsv);
title(tname);
exportgraphics(f,fname, 'resolution', 300);

end


