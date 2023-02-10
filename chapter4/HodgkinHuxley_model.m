% Tutorial 4.1, this code builds the Hodgkin-Huxley model as an oscillator
% By Clark Xu, Feb. 10, 2023

%% parameter
G_leak = 30e-9;
Gmax_Na = 12e-6;
Gmax_K = 3.6e-6;
E_Na = 0.045;
E_K = -0.082;
E_l = -0.060;
C_m = 100e-12;


%% time vector
tmax = 0.35;
dt = 1e-4;  % 0.1ms
tvec = 0:dt:tmax;
Nt = length(tvec);
Ntrial = 1;

%% vector for membrane potentials, m, h, and n
Vm = zeros(Ntrial, Nt);
mvec = zeros(size(Vm));
hvec = zeros(size(Vm));
nvec = zeros(size(Vm));

%% vector for Iapp
Iapp = ones(size(Vm));

Ibaseline = 0;
Iamplitude = 0.22e-9;

Iapp(1,:) = Iapp(1,:)*Ibaseline;
Iapp(1, 1:0.1/dt) = Iamplitude;



Vinit = -7.2e-3;
Vm(1,1) = Vinit;

for trial = 1:Ntrial
    
    for i = 2:Nt

        alphaM = (1e5 * (-Vm(trial,i-1)-0.045))/ ...
            (exp(100 * (-Vm(trial,i-1)-0.045))-1);
        betaM = 4e3 * exp((-Vm(trial,i-1)-0.070)/0.018);
        dm = alphaM*(1-mvec(trial,i-1)) - betaM*mvec(trial,i-1);
        mvec(trial,i) = mvec(trial,i-1)+dm*dt;

        alphaH = 70 * exp(50*(-Vm(trial,i-1)-0.070));    
        betaH = 1e3/(1+exp(100*(-Vm(trial,i-1)-0.040)));
        dh = alphaH*(1-hvec(trial,i-1)) - betaH*hvec(trial,i-1);
        hvec(trial,i) = hvec(trial,i-1)+dh*dt;

        alphaN = 1e4 * (-Vm(trial,i-1)-0.060)/ ...
            (exp(100*(-Vm(trial,i-1)-0.060))-1);
        betaN = 125 * exp((-Vm(trial,i-1)-0.070)/0.08);
        dn = alphaN*(1-nvec(trial,i-1)) - betaN*nvec(trial,i-1);
        nvec(trial,i) = nvec(trial,i-1)+dn*dt;

        dVm = G_leak * (E_l-Vm(trial,i-1)) + ...
            Gmax_Na * mvec(trial,i)^3 * hvec(trial,i) * (E_Na-Vm(trial,i-1)) + ...
            Gmax_K * nvec(trial,i)^4 * (E_K - Vm(trial,i-1)) + Iapp(trial,i);
        Vm(trial,i) = Vm(trial,i-1)+dVm*(dt/C_m);
    end
end

subplot(2,1,1);
plot(tvec,Iapp);

subplot(2,1,2);
plot(tvec, Vm);
