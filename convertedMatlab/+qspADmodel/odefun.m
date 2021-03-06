
function dcdt = odefun(t, c, p, de)
    % ODEFUN - system of equations of the AD model
    % Inputs:
    %   t - time 
    %   c - state variables
    %   p - parameter structure
    %   de - drug effects
    % Outputs:
    %   dcdt - Derivatives of the state variables
    % the model converted to matlab from original python script
    % by Jaimit Parikh
    
    
    %     # c : levels of 14 biological factors at t. ndarray of shape (14)
    SB = c(1); % skin barrier integrity
    IP = c(2); % infiltrated pathogens
    Th1 = c(3); % Th1
    Th2 = c(4); % Th2
    Th17 = c(5); % Th17
    Th22 = c(6); % Th22
    IL4 = c(7); % IL4
    IL13 = c(8); % IL13
    IL17 = c(9); % IL17
    IL22 = c(10); % IL22
    IL31 = c(11); % IL31
    IFNg = c(12); % IFNg
    TSLP = c(13); % TSLP
    OX40L = c(14); % OX40L
    %     # t : time (int)
    %     # de: effects of [placebo, IL4, IL13, IL17A, IL22, IL31, TSLP, OX40, rIFNg]
    %     # x : 51 parameter values   
   
    %     #   prepare parameter values
    k1 = p.k1; 
    k2 = p.k2; 
    k3 = p.k3;
    b1 = p.b1;
    b2 = p.b2;
    b3 = p.b3;
    b4 = p.b4;
    b5 = p.b5;
    d1 = p.d1;
    d2 = p.d2;
    d3 = p.d3;
    b6 = p.b6;
    d4 = p.d4;
    d5 = p.d5;
    d6 = p.d6;
    d7 = p.d7;
    b7 = p.b7;
    b8 = p.b8;
    d8 = p.d8;
    k5 = p.k5;
    k9 = p.k9;
    d9 = p.d9;
    b9 = p.b9;
    k6 = p.k6;
    k10 = p.k10;
    k7 = p.k7;
    k8 = p.k8;
    k11 = p.k11;
    k12 = p.k12;
    d10 = p.d10;
    k13 = p.k13;
    k14 = p.k14;
    d11 = p.d11;
    k15 = p.k15;
    k16 = p.k16;
    d12 = p.d12;
    k17 = p.k17;
    k18 = p.k18;
    d13 = p.d13;
    k19 = p.k19;
    k20 = p.k20;
    d14 = p.d14;
    k21 = p.k21;
    k22 = p.k22;
    d15 = p.d15;
    k23 = p.k23;
    k24 = p.k24;
    d16 = p.d16;
    k25 = p.k25;
    k26 = p.k26;
    d17 = p.d17;
    k3 = min(k3, de(1));
    ea2 = max(0.4396, de(10));
    k4 = d8;
    
    % #   effective concentration of cytokines (drug effects on cytokines)
    IL4e  = (1 - de(2))*IL4;
    IL13e = (1 - de(3)*ea2)*IL13;
    IL17e = (1 - de(4))*IL17;
    IL22e = (1 - de(5))*IL22;
    IL31e = (1 - de(6))*IL31;
    TSLPe = (1 - de(7))*TSLP;
    OX40Le = (1 - de(8))*OX40L;
    IFNge = IFNg + de(9);
    
    %#   ODEs
        
    %# skin barrier integrity        
    dSBdt = (1 - SB)*(k1 + k2*IL22e +k3) / ...
        ((1 + b1*IL4e)*(1 + b2*IL13e)*(1 + b3*IL17e)*(1 + b4*IL22e)*(1 + b5*IL31e))...
        - SB*(d1*(1 + d3*IP) + d2*IL31e);

    % # infiltrated pathogens        
    dIPdt = k4/(1 + b6*SB) -...
        IP*(((1 + d4*IP)*(1 + d5*IL17e)*(1 + d6*IL22e)*(1 + d7*IFNge))/...
        ((1 + b7*IL4e)*(1 + b8*IL13e)) + d8);

    % # Th cells
    dTh1dt = k5*IP*(1 + k9*IFNge) /(4 + k9*IFNge + k10*IL4e) - d9*Th1/(1 + b9*OX40Le); 
    dTh2dt = k6*IP*(1 + k10*IL4e) /(4 + k9*IFNge + k10*IL4e) - d9*Th2/(1 + b9*OX40Le);
    dTh17dt = k7*IP                    /(4 + k9*IFNge + k10*IL4e) - d9*Th17/(1 + b9*OX40Le);
    dTh22dt = k8*IP                    /(4 + k9*IFNge + k10*IL4e) - d9*Th22/(1 + b9*OX40Le);

     % #cytokines        
    dIL4dt  = k11*Th2 + k12 - d10*IL4;  %# IL4
    dIL13dt  = k13*Th2 + k14 - d11*IL13;  %# IL13
    dIL17dt  = k15*Th17 + k16 - d12*IL17;  %# IL17
    dIL22dt  = k17*Th22 + k18 - d13*IL22;  %# IL22
    dIL31dt = k19*Th2 + k20 - d14*IL31; %# IL31
    dIFNgdt = k21*Th1 + k22 - d15*IFNg; %# IFNg
    dTSLPdt = k23*IP + k24 - d16*TSLP; %# TSLP
    dOX40Ldt = k25*TSLPe + k26 - d17*OX40L; %# OX40L

    dcdt = [dSBdt, dIPdt, dTh1dt, dTh2dt, dTh17dt, dTh22dt,...
        dIL4dt, dIL13dt, dIL17dt, dIL22dt, dIL31dt, dIFNgdt, dTSLPdt, dOX40Ldt]';
end