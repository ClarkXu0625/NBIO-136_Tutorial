%% tutorial 4.3
% A two-compartment model of an intrinsically bursting neuron. 
% Written by Clark Xu, Feb. 17, 2023

%% parameters for simulation
% parameters are from table 4.7 and 4.8
A_S = 1/3;  % fraction of soma
A_D = 1-A_S;    % fraction of dendrite
Gs_leak = A_S*5e-9;     % somatic leak conductance
Gd_leak = A_D*5e-9;     % drendritic leak conductance
Gmax_Na = A_S*3e-6;     % maximum sodium conductance
Gmax_K = A_S*2e-6;  % maximum delayed rectifier conductance
Gmax_Ca = A_D*2e-6;     % maximum calcium conductance
Gmax_KCa = A_D*2.5e-6;  % max calcium-dependent potassium condunctance
Gmax_KAHP = A_D*40e-9;   % max after-hyperpolarization condunctance
G_link = 10e-9 ;    % link condunctance
E_Na = 0.060;   % sodium reveral potential
E_Ca = 0.080;   % calcium reversal potential
E_K = -0.075;   % potassium reversal potential
E_l = -0.060;   % leak reversal potential
C_S = A_S*100e-12;  % capacitance of soma
C_D = A_D*100e-12;  % capacitance of dendrite
Iapp_S = 0e-12;     % somatic applied current
Iapp_D = 0e-12;     % dendritic applied current
tao_Ca = 0.050;     % calcium decay time constant, 50ms
k = 5e6/A_D;    % conversion from charge to concentration


tmax = 2;
dt = 2e-6;  % timestep is 2us
tvec = 0:dt:tmax;
Nt = length(tvec);
Ntrial = 1;

%% create vectors for membrane potentials, gating variables, and Ca conc.
Vm_D = zeros(Ntrial, Nt);
Vm_S = zeros(size(Vm_D));
mvec = zeros(size(Vm_D));
nvec = zeros(size(Vm_D));
hvec = zeros(size(Vm_D));
mvec_Ca = zeros(size(Vm_D));
mvec_KCa = zeros(size(Vm_D));
mvec_KAHP = zeros(size(Vm_D));
Ca_conc_vec = zeros(size(Vm_D));

spk_detect=0;
spk_count=0;
up_th = -0.010;
down_th = -0.030;
%% start simulation
for i = 2:Nt
    [alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] ...
        = PR_soma_gating(Vm_S(1,i-1));
    [alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] ...
        = PR_dend_gating(Vm_D(1,i-1), Ca_conc_vec(1,i-1));

    % updating gating variables
    dm = gating_variable(alpha_m, beta_m, mvec(1,i-1), dt);
    mvec(1,i) = mvec(1,i-1)+dm;

    dn = gating_variable(alpha_n, beta_n, nvec(1,i-1), dt);
    nvec(1,i) = nvec(1,i-1)+dn;

    dh = gating_variable(alpha_h, beta_h, hvec(1,i-1), dt);
    hvec(1,i) = hvec(1,i-1)+dh;

    dm_Ca = gating_variable(alpha_mca, beta_mca, mvec_Ca(1,i-1), dt);
    mvec_Ca(1,i) = mvec_Ca(1,i-1)+dm_Ca;

    dm_KAHP = gating_variable(alpha_kahp, beta_kahp, mvec_KAHP(1,i-1), dt);
    mvec_KAHP(1,i) = mvec_KAHP(1,i-1)+dm_KAHP;

    dm_KCa = gating_variable(alpha_kca, beta_kca, mvec_KCa(1,i-1), dt);
    mvec_KCa(1,i) = mvec_KCa(1,i-1)+dm_KCa;
        
    % update Cacium concentration
    I_Ca = Gmax_Ca * (mvec_Ca(1,i))^2 * (E_Ca-Vm_D(1,i-1));
    dCa = (-Ca_conc_vec(1,i-1)/tao_Ca + k*I_Ca)*dt;
    Ca_conc_vec(1,i) = Ca_conc_vec(1,i-1)+dCa;
    X = min(4000*Ca_conc_vec(1,i),1);

    % update membrane potential for dendrite and soma
    dV_S = (Gs_leak * (E_l-Vm_S(1,i-1)) + ...
        Gmax_Na * (mvec(1,i)^2) * hvec(1,i) * (E_Na-Vm_S(1,i-1)) + ...
        Gmax_K * (nvec(1,i)^2) * (E_K-Vm_S(1,i-1)) + ...
        G_link * (Vm_D(1,i-1) - Vm_S(1,i-1)) + Iapp_S)*dt/C_S;
    Vm_S(1,i) = Vm_S(1,i-1)+dV_S;

    dV_D = (Gd_leak * (E_l-Vm_D(1,i-1)) + ...
        Gmax_Ca * (mvec_Ca(1,i))^2 * (E_Ca-Vm_D(1,i-1)) + ...
        Gmax_KCa * mvec_KCa(1,i) * X * (E_K-Vm_D(1,i-1)) + ...
        Gmax_KAHP * mvec_KAHP(1,i) * (E_K-Vm_D(1,i-1)) - ...
        G_link * (Vm_D(1,i-1) - Vm_S(1,i-1)) + Iapp_D)*dt/C_D;
    Vm_D(1,i) = Vm_D(1,i-1)+dV_D;

    if spk_detect
        if Vm_S(1,i)>up_th
            spk_detect=0;
            spk_count=spk_count+1;
        end
    else
        if Vm_S(1,i)<down_th
            spk_detect=1;
        end
    end

end

pks = findpeaks(Vm_S, 'MinPeakHeight', up_th, 'MinPeakProminence',up_th-down_th);
disp(pks)
disp(size(pks))
disp(spk_count)

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
plot(tvec,Vm_D);
figure(2)
plot(tvec,Vm_S);


function d = gating_variable(alpha, beta, previous, dt)
    eq = alpha/(alpha+beta);    % steady state
    tao = 1/(alpha+beta);   % time constant
    d = (eq-previous)*dt/tao;
end