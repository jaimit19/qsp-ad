% SENSPRCCREG - Example file by Jaimit Parikh to run the PRCC and Linear
% Regression based sensitivity analysis on the QSP AD model
% Need to be further generalized to run it on desired set of variables

import qspADmodel.*; % Import the functions from the qspADmodel package
import sa.*; % Import the functions from the sa (sensitivity analysis) package


% 1. Get parameters of the AD model
parameters = adModelPars();

% 2. Sample using Latin hypercube sampling the parameter values and running
% simulations for the samples

% LHS sampling for the desired parameters
parsName = ["k1",  "k2", "b1", "d1", "k5", "k6"];
nParameters = length(parsName);
problem.lowerBound = [2.6,  0.24,  6.3e-5,  0.15, 10, 25];
problem.upperBound = [10.4, 0.96, 2.52e-4, 0.6,  40, 100];
problem.num_vars = nParameters;
problem.names = parsName;
N = 1000; % number of samples
samples = latinHypercubeSampling(problem,  N);

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


% 3. Perform PRCC based sensitivity analysis
S_SB = prccAnalysis(samples, SB', "partialcorr");
S_Th1 = prccAnalysis(samples, Th1', "partialcorr");
plotSAresults(S_SB.prcc(1:end-1, end), S_Th1.prcc(1:end-1, end),...
    parsName, './figures/prcc.png', 'PRCC Analsis')


% 4. Perform Regression based sensitivity analysis

S_SB_reg = regAnalysis(samples, SB');
S_Th1_reg = regAnalysis(samples, Th1');
plotSAresults(S_SB_reg.S1, S_Th1_reg.S1, ...
    parsName, './figures/regAnalysis.png', 'Regression Analysis')


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


