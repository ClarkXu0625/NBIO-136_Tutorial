% Tutorial 3.2: Statistical properties of simulated spike trains. 
% Written by Clark Xu, 2/28/2023

% parameters
E_l = -0.070;
V_th = -0.050;
V_reset = -0.080;
del_th = 0.002;
G_l = 10e-9;
C_m = 100e-12;
a = 2e-9;
b = 0e-9;
tao_SRA = 0.150;

% time vector
dt = 1e-5;
tmax = 0.1;
tvec = 0:dt:tmax;
