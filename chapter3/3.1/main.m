% Tutorial 3.1, Generating receptive fields with spike-triggered averages
% Written by Clark Xu, last modified on 3/1/2023

%% parameters
E_l = -0.060;
V_th = -0.050;
V_reset = -0.080;
del_th = 0.002;
G_l = 8e-9;
C_m = 100e-12;
a = 10e-9;
b = 0.5e-9;
tau_SRA = 0.050;

%% set up vectors
Nblock = 4e5;       % number of time blocks
block = 5e-3;       % each time block is 5ms
Iapp_rand = (1e-9).*rand(1,Nblock) - 0.5e-9;  % applied current, uniformly random within ±0.5nA  
dt = 2e-5;          % 0.02ms, time-bin of neuron
new_dt = 0.001;     % 1ms, this analysis focus on changes on a time-scale of 1ms or more
tvec = 0:dt:(block*Nblock-dt);  % time vector for neuron simulation
Nt = length(tvec);              % Number of time steps for neuron simulation
Iapp = zeros(size(tvec));       % Applied current, same size as time vector
ratio = round(block/dt);        % times new vector greater than old

% 250 repetition for each randomly generized Iapp value
for i=1:Nblock
    Iapp(1, (i-1)*ratio+1 : i*ratio) = Iapp_rand(1,i);
end

% Set up vectors for membrane potential V for simulation, I_SRA, and spike 
% which records the time that the neuron fires
V = zeros(1, Nt);       % Membrane Potential
V(:,1)=E_l;             % initialize membrane potential
I_SRA = zeros(1,Nt);   % spike-rate adaptation current
spike = zeros(size(V));     % vector recording spike time

%% simulation starts here.
for i = 2:Nt-1
    dI_SRA = (a*(V(1,i-1)-E_l)-I_SRA(1,i-1))*dt/tau_SRA;
    I_SRA(1,i) = I_SRA(1,i-1)+dI_SRA;

    dV = (G_l*(E_l-V(1,i-1)+del_th*exp((V(1,i-1)-V_th)/del_th))- ...
        I_SRA(1,i)+Iapp(1,i))*dt/C_m;
    V(1,i) = V(1,i-1)+dV;
            
    % spike occurs
    if V(1,i)>V_th
        V(1,i)=V_reset;  % reset membrane potential
        spike(1, i) = 1;
        % increment I_SRA
        I_SRA(1,i)=I_SRA(1,i)+b;
                        
    end    
end

Iapp = expandbin(Iapp, dt, new_dt);     % expand Iapp timebin
spike = expandbin(spike, dt, new_dt);   % expand spike timebin
spike(find(spike)) = 1;   % replace fractions with a "1" for a spike
[sta, tcorr] = STA(Iapp, spike, new_dt);   


figure(1)
hold on
plot(tcorr,sta)