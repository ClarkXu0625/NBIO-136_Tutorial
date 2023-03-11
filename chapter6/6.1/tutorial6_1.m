% Tutorial 6.1. Bistability and oscillations in a firing-rate model with 
% feedback. 
% Written by Clark Xu, last modified on 3/10/2023

question_number = 1;
question = 'a';

%% parameter
x = 1.2;
t_0 = 0.1;
r_max = 100;
sigma = 0.5;
tau_r = 0.010;  % time constant for changes in firing rate


%% vectors
tmax = 20;  % 20 second simulation
dt = 1e-5;  % 0.01ms
tvec = 0:dt:tmax;
Nt = length(Nt);    % number of time steps
r = zeros(1,Nt);
D = zeros(size(r));  % depression variable
s = zeros(size(r));  % synaptic gating variable
s_in = zeros(size(r));

switch question_number
    case 1
        alpha_0 = 0.5;
        W_EE = 8;
        p_r = 1;
        tau_s = 0.002;
        s_0 = 0;    % initial condition for s
        r_0 = 0;    % initial condition for r
        s_in(10/dt:(10+0.050)/dt) = 0.05;   % set temporary input strength
        D(:) = 1;
end

%% start simulation
for i = 2:Nt
    s(i) = s(i-1) + (-s(i-1)/tau_s+alpha_0*D(i-1)*p_r*r(i-1*(1-s(i-1))))*dt;
    
    S = W_EE*s(i)+s_in(i);
    
    f = firing_rate_curve(r_0, r_max, sigma, S);

    r(i) = r(i-1) + (-r(i-1)+f)/tau_r*dt;
        rem
end

s_value = 0:0.01:1;
f = firing_rate_curve(r_0, r_max, sigma, s*W_EE);
figure(1)
plot(s_value, f);


