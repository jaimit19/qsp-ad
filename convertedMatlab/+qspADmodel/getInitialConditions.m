function IC = getInitialConditions()
% Returns array of initial conditions of the AD model

s_0 = 0.5931;
p_0 = 0.4069;
Th1_0 = 3.1;
Th2_0 = 8.7;
Th17_0 = 2.0;
Th22_0 = 21.0;
IL4_0 = 38.0;
IL13_0 = 40.5;
IL17_0 = 5.4;
IL22_0 = 3.0;
IL31_0 = 2.0;
IFNg_0 = 1.5;
TSLP_0 = 4.4;
OX40_0 = 9.7;
    
IC = [s_0, p_0, Th1_0, Th2_0, Th17_0, Th22_0, IL4_0, ...
    IL13_0, IL17_0, IL22_0, IL31_0, IFNg_0, TSLP_0, OX40_0];
end