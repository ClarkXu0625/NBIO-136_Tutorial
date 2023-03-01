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
tau_SRA = 0.150;
sigma = 50e-12;


% time vector
dt = 1e-5;
tmax = 0.1;
tvec = 0:dt:tmax;
Nt = length(tvec);
Iapp = normrnd(0,sigma/sqrt(dt),1,Nt);

% Set up vectors for membrane potential V for simulation, I_SRA, and spike 
% which records the time that the neuron fires
V = zeros(1, Nt);       % Membrane Potential
V(:,1)=E_l;             % initialize membrane potential
I_SRA = zeros(size(V));   % spike-rate adaptation current
spikes = zeros(size(V));
noise_vec = randn(size(V))*sigma*sqrt(dt);

%% simulation starts here.
for i = 2:Nt-1
    
    dI_SRA = (a*(V(1,i-1)-E_l)-I_SRA(1,i-1))*dt/tau_SRA;
    I_SRA(1,i) = I_SRA(1,i-1)+dI_SRA;

    dV = (G_l*(E_l-V(1,i-1)+del_th*exp((V(1,i-1)-V_th)/del_th))- ...
        I_SRA(1,i)+Iapp(1,i))*dt/C_m;
    V(1,i) = V(1,i-1)+dV+noise_vec(i);
            
    % spike occurs
    if V(1,i)>V_th
        V(1,i)=V_reset;             % reset membrane potential      
        I_SRA(1,i)=I_SRA(1,i)+b;    % increment I_SRA  
        spikes(1,i) = 1;            % record spike time
    end    
end

spike_index = find(spikes);
ISI = zeros(1, length(spike_index)-1);
for i = 1:length(ISI)
    ISI(i) = spike_index(i+1)-spike_index(i);
end

plot(tvec, Iapp)
disp(spike_index)

histogram(ISI,25)
cv = std(ISI)/mean(ISI);
disp(cv)
