% Tutorial 3.2: Statistical properties of simulated spike trains. 
% Written by Clark Xu, last modified in 3/3/2023

clear

%% parameters
E_l = -0.070;
V_th = -0.050;
V_reset = -0.080;
del_th = 0.002;
G_l = 10e-9;
C_m = 100e-12;
a = 2e-9;
b = 1e-9;
tau_SRA = 0.150;
sigma = 50e-12;
Iapp_cons = 0e-9;

% time vector
dt = 1e-5;  % 0.01ms time step
tmax = 100;  % duration of 100 seconds
tvec = 0:dt:tmax;
Nt = length(tvec);

% Iapp has mean of 0, but with standard deviation of sigma/sqrt(dt)
%Iapp = normrnd(0,sigma/sqrt(dt),[1,Nt]);  
Iapp = randn(size(tvec))*sigma/sqrt(dt) + ones(size(tvec))*Iapp_cons;


% Set up vectors for membrane potential V for simulation, I_SRA, and spike 
% which records the time that the neuron fires
V = zeros(1, Nt);       % Membrane Potential
V(:,1)=E_l;             % initialize membrane potential
I_SRA = zeros(size(V));   % spike-rate adaptation current
spikes = zeros(size(V));    % vector records when does neuron fires
%noise_vec = randn(size(V))*sigma*sqrt(dt);

%% simulation starts here.
for i = 2:Nt-1

    % spike occurs
    if V(1,i)>V_th
        V(1,i)=V_reset;             % reset membrane potential      
        I_SRA(1,i)=I_SRA(1,i)+b;    % increment I_SRA  
        spikes(1,i) = 1;            % record spike time
    end   

    dI_SRA = (a*(V(1,i)-E_l)-I_SRA(1,i))*dt/tau_SRA;
    I_SRA(1,i+1) = I_SRA(1,i)+dI_SRA;

    dV = (G_l*(E_l-V(1,i)+del_th*exp((V(1,i)-V_th)/del_th))- ...
        I_SRA(1,i)+Iapp(1,i))*dt/C_m;
    V(1,i+1) = V(1,i)+dV;      
end

spike_index = find(spikes);     % list of index when neruon fires
ISI = zeros(1, length(spike_index)-1);
min =10;

% ISI records the ISI for each spike, in unit of second
for i = 1:length(ISI)
    ISI(i) = (spike_index(i+1)-spike_index(i))*dt;
    if ISI(i)<min
        min=ISI(i);
        time=i;
    end
end
disp(min)
disp(min/dt)

figure(1)
histogram(ISI*1000,[5:25:355])


xlabel("ISI (ms)")
ylabel("Count")


bin_size = 0.010: 0.010: 1;
[fano, variance, mean] = fano_factor(spikes, dt, bin_size);
figure(2)
subplot(3,1,1)
plot(bin_size,fano);
xlabel("bin size (s)")
ylabel("fano factor")

subplot(3,1,2)
plot(bin_size, mean);

subplot(3,1,3)
plot(bin_size, variance)

figure(10)
subplot(3,1,1)
plot(tvec, Iapp)
subplot(3,1,2)
plot(tvec, V*1000)
ylabel("membrane potential")
subplot(3,1,3)
plot(tvec,spikes)
axis([0 tmax 0 2]);
