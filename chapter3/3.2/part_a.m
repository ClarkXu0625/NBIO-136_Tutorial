% Tutorial 3.2: Statistical properties of simulated spike trains. 
% Written by Clark Xu, last modified in 3/9/2023

clear

question_number = 'c';
disp("question " + question_number+")")
%% parameters
E_l = -0.070;
V_th = -0.050;
V_reset = -0.080;
del_th = 0.002;
G_l = 10e-9;
C_m = 100e-12;
a = 2e-9;
tau_SRA = 0.150;

switch question_number
    case 'a'
        b = 0e-9;
        sigma = 50e-12;
        Iapp_cons = 0e-9;
    case 'b'
        b = 1e-9;
        sigma = 50e-12;
        Iapp_cons = 0e-9;
    case {'c'}
        b = 0e-9;
        sigma = 20e-12;
        Iapp_cons = 0.1e-9;
end

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

spike_index = dt*find(spikes);     % list of index when neruon fires (unit = second)
ISI=1000*diff(spike_index);             % ISIs converted to ms
cv = std(ISI,1)/mean(ISI);      % CV for inter-spike interval

% calculate cv2
sum2=0.0;
for i=2:length(ISI)
    sum2 = sum2 + sqrt((ISI(i)-ISI(i-1))*(ISI(i)-ISI(i-1))) ...
        /(ISI(i)+ISI(i-1))*(2.0);
end
cv2 = sum2 /(length(ISI)-1);

min_ISI = min(ISI);
max_ISI = max(ISI);
%disp(min_ISI)


%% plot
figure(1)
histogram(ISI, min_ISI : (max_ISI-min_ISI)/25 : max_ISI); % 25 time bins
xlabel("ISI (ms)")
ylabel("Count")
title("cv = " + num2str(cv) + " cv2 = " + num2str(cv2))

switch question_number 
    case {a, b}
        %% find the fano factor by bin size, varying from 10ms to 1s.
        bin_size = 0.010: 0.001: 1;
        [fano, ~, ~] = fano_factor(spikes, dt, bin_size);
        figure(2)
        plot(bin_size,fano);
        xlabel("bin size (s)")
        ylabel("fano factor")
end

[fano, variance, mean] = fano_factor(spikes, dt, 0.1);
disp("  When time-bin is 100ms: " )
disp("    Variance="+num2str(variance))
disp("    Mean=" + num2str(mean))
disp("    Fano_factor="+ num2str(fano))


