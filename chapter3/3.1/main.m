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
Nblock = 4e5;   % number of time blocks
block = 5e-3;   % each time block is 5ms

% applied current, random within ±0.5nA
Iapp_rand = (1e-9).*rand(1,Nblock) - 0.5e-9; 
dt = 2e-5;  % 0.02ms, time-bin of neuron
new_dt = 0.001; % 1ms, only interested in changes on a time-scale of 1ms or more
tvec = 0:dt:(block*Nblock-dt); 
Nt = length(tvec);
Iapp = zeros(size(tvec));
ratio = round(block/dt);   % times new vector greater than old

for i=1:Nblock
    Iapp(1, (i-1)*ratio+1 : i*ratio) = Iapp_rand(1,i);
end


Ntrial = 1;
% Set up vectors for Vm, G_SRA, and I_SRA
V = zeros(Ntrial, Nt);  % Membrane Potential
V(:,1)=E_l; % initialize membrane potential
I_SRA = zeros(Ntrial,Nt);
spike = zeros(size(V));

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

Iapp = expandbin(Iapp, dt, new_dt);
spike = expandbin(spike, dt, new_dt);
spike(find(spike)) = 1;   % replace fractions with a "1" for a spike
[sta, tcorr] = STA(Iapp, spike, new_dt);

figure(2)
plot(tcorr,sta)