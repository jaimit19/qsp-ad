% RUNADMODEL - Example file by Jaimit Parikh to run the AD QSP model  
% https://pubmed.ncbi.nlm.nih.gov/33894014/converted to matlab
% from the original python scripts https://github.com/Tanaka-Group/AD_QSP_model
% 

import qspADmodel.*; % Import the functions from the qspADmodel package

parameters = adModelPars(); % Load parameters of the AD model
de = adDrugEffects(); % Set and load the drug effects to be simulated
IC = getInitialConditions(); % get the initial values of the state variables
t0 = 0; tfinal = 1000; % set the time for simulation

% Simulate the model 
[T,Y] = ode15s(@(t, y)odefun(t, y, parameters, de),...
    [t0 tfinal], IC);

% Plot the state variables with time and compare with original python
% results
plotAD(T, Y, './figures/stateVariables.png');


function plotAD(T, Y, fname)
f = figure('DefaultAxesFontSize',12,...
    'Position', [20 20 1800 900]);
sp = ["SB", "IP", "Th1", "Th2", "Th17", "Th22", "IL4", ...
    "IL13", "IL17", "IL22", "IL31", "IFNg", "TSLP", "OX40L"];
tb =  readtable('python_res.csv');
for idx = 1:size(Y, 2)
    subplot(3, 5, idx);
    plot(T, Y(:, idx), 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'matlab');
    hold on ;
    plot(linspace(0, 1000, height(tb)), tb(:, idx).Variables, 'Color', 'r', 'LineWidth', 2, ...
        'LineStyle', '--', 'DisplayName', 'python');
    hold off;
    xlabel('Time, weeks');
    ylabel(sp(idx));
    xlim([0, 50]);
end
legend();
exportgraphics(f, fname, 'resolution', 300);
end

function e = EASI(sim)
    s = sim(:,0);
    p = sim(:,1);
    e = 72 * (2*p + 2*(1-s)) / 4;
end
