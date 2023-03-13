% Tutorial 6.1. Bistability and oscillations in a firing-rate model with 
% feedback. 
% Written by Clark Xu, last modified on 3/10/2023
clear

question_number = 2;
question = 'a';

%% parameter
x = 1.2;    % exponent
r_0 = 0.1;   % the rate with no inputs
r_max = 100;
sigma = 0.5;
tau_r = 0.010;  % time constant for changes in firing rate
start_time = 10;    % start time of temp input strength

%% vectors
tmax = 20;  % 20 second simulation
dt = 1e-4;  % 0.1ms
tvec = 0:dt:tmax;
Nt = length(tvec);    % number of time steps
r = zeros(1, Nt);
D = zeros(size(r));  % depression variable
s = zeros(size(r));  % synaptic gating variable
s_in = zeros(size(r));


switch question_number
    case 1
        alpha_0 = 0.5;
        W_EE = 8;
        p_r = 1;
        tau_s = 0.002;
        duration = 0.050;
        amplitude = 0.05;
        D = ones(size(r));
    case 2
        alpha_0 = 0.5;
        W_EE = 60;
        p_r = 0.2;
        tau_s = 0.002;
        tau_D = 0.250;
        duration = 2;
        amplitude = 0.002;
    case 3
        alpha_0 = 0.5;
        p_r = 0.5;
        W_EE = 35;
        r_0 = -0.1;
        tau_s = 0.002;
        tau_D = 0.250;
        duration = 0.050;
        amplitude = 0.05;
        if question == 'c'
            r(1) = 9;
            D(1) = 1/(1+p_r*r(1)*tau_D);
            s(1) = alpha_0*D(1)*p_r*r(1)*tau_s/(1+alpha_0*D(1)*p_r*r(1)*tau_s);
        end
    case 4
        alpha_0 = 0.25;
        p_r = 1;
        W_EE = 35;
        r_0 = -0.1;
        tau_s = 0.002;
        tau_D = 0.125;
        duration = 0.050;
        amplitude = 0.05;
end

% set temporary input strength
s_in(round(start_time/dt):round((start_time+duration)/dt)) = amplitude;   

%% start simulation
for i = 2:Nt
    % update depression variable if D is not fixed
    if (exist('tau_D'))
        D(i) = D(i-1) + ((1-D(i-1))/tau_D - p_r*D(i-1)*r(i-1))*dt;
    end

    % update synaptic gating variable
    s(i) = s(i-1) + (-s(i-1)/tau_s+alpha_0*D(i-1)*p_r*r(i-1)*(1-s(i-1)))*dt;
    
    % total synaptic input
    S = W_EE*s(i)+s_in(i);
    
    % get firing rate
    f = firing_rate_curve(r_0, r_max, x, sigma, S);

    % firing rate
    r(i) = r(i-1) + (-r(i-1)+f)/tau_r*dt;

    % 
    r(i) = min(r(i), r_0+r_max);
    r(i) = max(r(i), 0);
    
end

%% for question (a), plot f_curve and s_curve
s_value = 0:0.01:1;
f_curve = firing_rate_curve(r_0, r_max, x, sigma, s_value*W_EE);
figure(1)
subplot(2,1,1)
plot(s_value,f_curve)

s_curve = gating_variable_s(0:1:r_max, p_r, alpha_0, tau_s);
subplot(2,1,2)
plot(0:1:r_max,s_curve)


figure(2)
plot(tvec, r);
