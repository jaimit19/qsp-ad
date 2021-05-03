# -*- coding: utf-8 -*-
"""
This software is released under the MIT License, see LICENSE.txt.
Copyright (c) 2021 Takuya Miyano

QSP model simulates %improved EASI and EASI-75

"""

import numpy as np
from scipy.integrate import odeint
from numba import jit
from matplotlib import pyplot as plt
import pandas as pd

class ODE(object):
    def __init__(self, diff_eq, init_con):
        self.diff_eq  = diff_eq
        self.init_con = init_con
            
    def cal_equation(self, t_end, drug_effect, x):
        dt = 0.1 # delta time (week)
        N = round(t_end/dt) + 1 #Time steps
        t = np.linspace(0, t_end, N) # prepare time
        v = odeint(self.diff_eq, self.init_con, t, rtol=1e-8, atol=1e-6, args=(drug_effect, x))
        return v
    
@jit('f8[:](f8[:,:])', nopython=True)
def EASI(sim):
    s = sim[:,0]
    p = sim[:,1]
    e = 72 * (2*p + 2*(1-s)) /4
    return e
    
@jit('f8[:](f8[:],f8,f8[:],f8[:])', nopython=True, debug=True)
def diff_eq(c, t, de, x):
# c : levels of 14 biological factors at t. ndarray of shape (14)
#    c[0] : skin barrier integrity
#    c[1] : infiltrated pathogens
#    c[2] : Th1
#    c[3] : Th2
#    c[4] : Th17
#    c[5] : Th22
#    c[6] : IL4
#    c[7] : IL13
#    c[8] : IL17
#    c[9] : IL22
#    c[10]: IL31
#    c[11]: IFNg
#    c[12]: TSLP
#    c[13]: OX40L
# t : time (int)
# de: effects of [placebo, IL4, IL13, IL17A, IL22, IL31, TSLP, OX40, rIFNg]
# x : 51 parameter values   

#   prepare parameter values
    k1 = x[0]
    k2 = x[1]
    k3 = min([x[2],de[0]])
    b1 = x[3]
    b2 = x[4]
    b3 = x[5]
    b4 = x[6]
    b5 = x[7]
    d1 = x[8]
    d2 = x[9]
    d3 = x[10]
    b6 = x[11]
    d4 = x[12]
    d5 = x[13]
    d6 = x[14]
    d7 = x[15]
    b7 = x[16]
    b8 = x[17]
    d8 = x[18]
    k5 = x[19]
    k9 = x[20]
    d9 = x[21]
    b9 = x[22]
    k6  = x[23]
    k10 = x[24]
    k7  = x[25]
    k8 =  x[26]
    k11 = x[27]
    k12 = x[28]
    d10 = x[29]
    k13 = x[30]
    k14 = x[31]
    d11 = x[32]
    k15 = x[33]
    k16 = x[34]
    d12 = x[35]
    k17 = x[36]
    k18 = x[37]
    d13 = x[38]
    k19 = x[39]
    k20 = x[40]
    d14 = x[41]
    k21 = x[42]
    k22 = x[43]
    d15 = x[44]
    k23 = x[45]
    k24 = x[46]
    d16 = x[47]
    k25 = x[48]
    k26 = x[49]
    d17 = x[50]
    ea2 = max([0.4396, de[9]])
    k4 = d8

#   effective concentration of cytokines (drug effects on cytokines)
    IL4  = (1 - de[1])*c[6]
    IL13 = (1 - de[2]*ea2)*c[7]
    IL17 = (1 - de[3])*c[8]
    IL22 = (1 - de[4])*c[9]
    IL31 = (1 - de[5])*c[10]
    TSLP = (1 - de[6])*c[12]
    OX40 = (1 - de[7])*c[13]
    IFNg = c[11] + de[8]

#   ODEs
    # skin barrier integrity        
    dc0dt = (1 - c[0])*(k1 + k2*IL22 +k3)/((1 + b1*IL4)*(1 + b2*IL13)*(1 + b3*IL17)*(1 + b4*IL22)*(1 + b5*IL31)) - c[0]*(d1*(1 + d3*c[1]) + d2*IL31)

    # infiltrated pathogens        
    dc1dt = k4/(1 + b6*c[0]) - c[1]*(((1 + d4*c[1])*(1 + d5*IL17)*(1 + d6*IL22)*(1 + d7*IFNg))/((1 + b7*IL4)*(1 + b8*IL13)) + d8)

    # Th cells
    dc2dt = k5*c[1]*(1 + k9*IFNg) /(4 + k9*IFNg + k10*IL4) - d9*c[2]/(1 + b9*OX40) # Th1
    dc3dt = k6*c[1]*(1 + k10*IL4) /(4 + k9*IFNg + k10*IL4) - d9*c[3]/(1 + b9*OX40) # Th2
    dc4dt = k7*c[1]               /(4 + k9*IFNg + k10*IL4) - d9*c[4]/(1 + b9*OX40) # Th17
    dc5dt = k8*c[1]               /(4 + k9*IFNg + k10*IL4) - d9*c[5]/(1 + b9*OX40) # Th22

    # cytokines        
    dc6dt  = k11*c[3] + k12 - d10*c[6]  # IL4
    dc7dt  = k13*c[3] + k14 - d11*c[7]  # IL13
    dc8dt  = k15*c[4] + k16 - d12*c[8]  # IL17
    dc9dt  = k17*c[5] + k18 - d13*c[9]  # IL22
    dc10dt = k19*c[3] + k20 - d14*c[10] # IL31
    dc11dt = k21*c[2] + k22 - d15*c[11] # IFNg
    dc12dt = k23*c[1] + k24 - d16*c[12] # TSLP
    dc13dt = k25*TSLP + k26 - d17*c[13] # OX40L

    dcdt = np.array([dc0dt, dc1dt, dc2dt, dc3dt, dc4dt, dc5dt, 
                     dc6dt, dc7dt, dc8dt, dc9dt, dc10dt, dc11dt, dc12dt, dc13dt])
    return dcdt

    
def simulate_one(x):    
    # simulate steady state (1000 weeks) /baseline levels
    # initial conditions for simulating steady-state levels of biological factors 
    s_0    = np.float64(0.5931)
    p_0    = np.float64(0.4069)
    Th1_0  = np.float64(3.1)
    Th2_0  = np.float64(8.7)
    Th17_0 = np.float64(2.0)
    Th22_0 = np.float64(21.0)
    IL4_0  = np.float64(38.0)
    IL13_0 = np.float64(40.5)
    IL17_0 = np.float64(5.4)
    IL22_0 = np.float64(3.0)
    IL31_0 = np.float64(2.0)
    IFNg_0 = np.float64(1.5)
    TSLP_0 = np.float64(4.4)
    OX40_0 = np.float64(9.7)
    init_cond = np.array([s_0, p_0, Th1_0, Th2_0, Th17_0, Th22_0, IL4_0, 
                          IL13_0, IL17_0, IL22_0, IL31_0, IFNg_0, TSLP_0, OX40_0], dtype='float64')

    ode = ODE(diff_eq, init_cond)
    drug_effect = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float64')
    sim_0 = ode.cal_equation(1000, drug_effect, x) # 1000 days for steady state
    init_cond2 = sim_0[10000,:]

    # use steady-state level as baseline levels (initial condition)
    ode = ODE(diff_eq, init_cond2)
    T_end = 24 # weeks
    # Placebo (other)
    drug_effect = np.array([1E20, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float64')
    sim_1 = ode.cal_equation(T_end, drug_effect, x)
    b = np.concatenate([np.array([EASI(sim_0)[10000]]), sim_0[10000,:]])
    Res_series = np.concatenate([b, \
                              EASI(sim_1)])
    return [sim_0, sim_1], Res_series

if __name__ == "__main__":
    mu = np.loadtxt("mu.csv", delimiter = ",", dtype = float).reshape(-1,1)
    sigma = np.loadtxt("sigma.csv", delimiter = ",", dtype = float).reshape(-1,1)
    n_patients = 1

    random_list = np.random.randn(51, n_patients)
    virtual_subjects = random_list*np.abs(sigma) + mu
    virtual_subjects = np.exp(virtual_subjects)
    sampleList = [(i, virtual_subjects) for i in range(n_patients)]
    
    x = np.array([5.29667675e+00, 4.86649030e-01, 7.76965391e+00, 1.26193679e-04,
       7.43404696e-03, 3.08019761e-02, 8.52311950e-01, 2.93506777e-01,
       2.90507430e-01, 1.47857174e-01, 8.05256065e+00, 1.25236637e+00,
       1.32631235e-01, 2.69537776e-03, 6.68590530e-03, 6.59076102e-05,
       3.23512783e-04, 2.98628946e-02, 2.18016957e-01, 1.99142857e+01,
       2.68748146e-02, 5.32373474e-01, 1.28528412e-01, 4.90428220e+01,
       3.23308501e-03, 1.39889673e+01, 7.10889091e+01, 2.69677906e+02,
       1.21634465e+04, 4.93984515e+02, 1.73131251e+03, 8.95427238e+03,
       3.51273458e+02, 5.47583379e+01, 8.57494380e+00, 1.31916247e+01,
       3.18021649e+00, 4.47615740e+00, 4.54857072e+01, 4.02201776e+00,
       4.31419500e+00, 2.41647496e+01, 1.20454230e+00, 9.76753563e+00,
       1.37577244e+01, 4.79706180e+01, 2.94988412e+01, 2.06391084e+01,
       1.04394248e+00, 2.59172402e+00, 1.29874779e+00])
    
    [sim_0, sim_1], Res_series = simulate_one(x)
    pd.DataFrame(sim_0).to_csv('python_res.csv', header=False, index=False)