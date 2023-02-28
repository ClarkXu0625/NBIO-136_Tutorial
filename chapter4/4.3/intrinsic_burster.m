%% tutorial 4.3
% A two-compartment model of an intrinsically bursting neuron. 
% Written by Clark Xu. Last modify on 2/27/2023.

question_number = 2;    % allow 2, 3 (4), 5, 6
mode='d';   % switch to either 's' or 'd' in q6.
disp("question number "+num2str(question_number) + ", mode=" +mode)

% vector for different values of membrane potential
Vm_values = -0.085:0.001:0.050;

% vector of calcium concentration
Ca_conc = 0: 5e-5: 2e-3;

[alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] = PR_soma_gating(Vm_values);
[alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] = PR_dend_gating( Vm_values, Ca_conc );

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%% Question 2, 
% Plot all 12 of the rate constants for gating variables 
% Figure 10 includes the gating variables dependent on membrane potential 
% or Calcium concentration
switch question_number
    case 2
        figure(1)
        
        subplot(2,3,1)
        plot(Vm_values*1000,alpha_m)
        hold on
        plot(Vm_values*1000,beta_m)
        legend(["alpha_m", "beta_m"])
        ylabel("Rate constant for m")
        
        subplot(2,3,2)
        plot(Vm_values*1000,alpha_h)
        hold on
        plot(Vm_values*1000,beta_h)
        legend(["alpha_h", "beta_h"])
        ylabel("Rate constant for h")
        
        subplot(2,3,3)
        plot(Vm_values*1000,alpha_n)
        hold on
        plot(Vm_values*1000,beta_n)
        legend(["alpha_n", "beta_n"])
        ylabel("Rate constant for n")
        
        subplot(2,3,4)
        plot(Vm_values*1000, alpha_mca)
        hold on
        plot(Vm_values*1000, beta_mca)
        legend(["alpha_mca", "beta_mca"])
        ylabel("Rate constant for m_Ca")
        
        subplot(2,3,5)
        plot(Vm_values*1000, alpha_kca)
        hold on
        plot(Vm_values*1000, beta_kca)
        legend(["alpha_kca", "beta_kca"])
        ylabel("Rate constant for m_KCa")
        
        % add label 
        for i=1:5
            subplot(2,3,i)    
            xlabel("Membrane potential (mV)")
        end
        
        subplot(2,3,6)
        plot(Ca_conc, alpha_kahp)
        hold on
        plot(Ca_conc, beta_kahp)
        legend(["alpha_kahp", "beta_kahp"])
        ylabel("Rate constant for m_KAHP")
        xlabel("Calcium concentration (M)")
        
        sgtitle("Rate Constant Relationship with Membrane Potential (mV) or Ca" + ...
            "Concentration (M)", 'FontWeight', 'Bold')
end

%% parameters for simulation for question 3-6
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
G_link = 50e-9 ;    % link condunctance
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

% list of G_link and Iapp_D values, applied in question 5 and 6 to explore
% their relationship with spike count.
G_link_values = 0:5:100;
Iapp_values = 0:10:200;


%% Set up time vector
tmax = 2;
dt = 2e-6;  % timestep is 2us
tvec = 0:dt:tmax;
Nt = length(tvec);
Ntrial = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up multiple trials for qeustion 5 and 6, exploring relationships
switch question_number       
    case 5
        Ntrial = length(G_link_values)+1;
        spk = zeros(1,Ntrial-1);    % number of spikes
        burst = zeros(1, Ntrial-1); % number of bursts
        rate = zeros(1, Ntrial-1);  % average number of spikes per burst

        % the parameters that are used in plotting for question 5
        G_link = 15e-9; 
        Iapp_D = 0;
    case 6
        Ntrial = length(Iapp_values)+1;
        spk = zeros(1,Ntrial-1);
        burst = zeros(1,Ntrial-1);
        rate = zeros(1, Ntrial-1);

        % the parameters that are used in plotting for question 4
        G_link = 50e-9;
        Iapp_D = 100e-12;
        Iapp_S = 0e-12;
end

%% Create vectors for dynamic parameters, 
% compartment membrane potentials, gating variables, and Ca concentration
Vm_D = zeros(1, Nt);    % dendritic membrane potential
Vm_S = zeros(size(Vm_D));   % somatic membrane potential
mvec = zeros(size(Vm_D));   % gating variable m
nvec = zeros(size(Vm_D));   % gating variable n
hvec = zeros(size(Vm_D));   % gating variable h
mvec_Ca = zeros(size(Vm_D));    % gating variable m_Ca
mvec_KCa = zeros(size(Vm_D));   % gating variable m_KCa
mvec_KAHP = zeros(size(Vm_D));  % gating variable m_KAHP
Ca_conc_vec = zeros(size(Vm_D));    % concentration


up_th = -0.010;     % membrane potential exceeding up_threshold counts as a spike when spk_detect avaiable
down_th = -0.030;   % when membrane potential drops below down_threhold, spk_detect becomes available
min_interval = 0.012;

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');


%% start simulation
for trial = 1:Ntrial
    % initialize spk_detect and spk_count for each trial
    spk_detect=0;   % check if spike is available to detect
    spk_count=0;    % count number of spikes during the simulation
    burst_count = 0;    % number of burst in this trial
    last_spike = 0; % time for last spike

    % simulate different parameters values from the second trial.
    if question_number==5 && (trial>1)
        G_link = G_link_values(trial-1)*1e-9;
    elseif question_number==6 && (trial>1)
        if mode=='d'
            Iapp_D = Iapp_values(trial-1)*1e-12;
            Iapp_S = 0;
        elseif mode=='s'
            Iapp_S = Iapp_values(trial-1)*1e-12;
            Iapp_D = 0;
        end
    end

    for i = 2:Nt
    
        [alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] ...
            = PR_soma_gating(Vm_S(1,i-1));
        [alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] ...
            = PR_dend_gating(Vm_D(1,i-1), Ca_conc_vec(1,i-1));
    
        %% update gating variables
        % gating variable m
        dm = gating_variable(alpha_m, beta_m, mvec(1,i-1), dt);
        mvec(1,i) = mvec(1,i-1)+dm;
    
        % gating variable n
        dn = gating_variable(alpha_n, beta_n, nvec(1,i-1), dt);
        nvec(1,i) = nvec(1,i-1)+dn;
    
        % gating variable h
        dh = gating_variable(alpha_h, beta_h, hvec(1,i-1), dt);
        hvec(1,i) = hvec(1,i-1)+dh;
    
        % gating variable m_Ca
        dm_Ca = gating_variable(alpha_mca, beta_mca, mvec_Ca(1,i-1), dt);
        mvec_Ca(1,i) = mvec_Ca(1,i-1)+dm_Ca;
    
        % gating variable m_KAHP
        dm_KAHP = gating_variable(alpha_kahp, beta_kahp, mvec_KAHP(1,i-1), dt);
        mvec_KAHP(1,i) = mvec_KAHP(1,i-1)+dm_KAHP;
    
        % gating variable m_KCa
        dm_KCa = gating_variable(alpha_kca, beta_kca, mvec_KCa(1,i-1), dt);
        mvec_KCa(1,i) = mvec_KCa(1,i-1)+dm_KCa;
            
        %% update Cacium concentration
        I_Ca = Gmax_Ca * (mvec_Ca(1,i))^2 * (E_Ca-Vm_D(1,i-1)); % Ca current
        dCa = (-Ca_conc_vec(1,i-1)/tao_Ca + k*I_Ca)*dt;
        Ca_conc_vec(1,i) = Ca_conc_vec(1,i-1)+dCa;
        X = min(4000*Ca_conc_vec(1,i),1);
    
        %% update membrane potentials
        % update somatic membrane potential
        dV_S = (Gs_leak * (E_l-Vm_S(1,i-1)) + ...
            Gmax_Na * (mvec(1,i)^2) * hvec(1,i) * (E_Na-Vm_S(1,i-1)) + ...
            Gmax_K * (nvec(1,i)^2) * (E_K-Vm_S(1,i-1)) + ...
            G_link * (Vm_D(1,i-1) - Vm_S(1,i-1)) + Iapp_S)*dt/C_S;
        Vm_S(1,i) = Vm_S(1,i-1)+dV_S;
    
        % update dendritic membrane potential
        dV_D = (Gd_leak * (E_l-Vm_D(1,i-1)) + ...
            Gmax_Ca * (mvec_Ca(1,i))^2 * (E_Ca-Vm_D(1,i-1)) + ...
            Gmax_KCa * mvec_KCa(1,i) * X * (E_K-Vm_D(1,i-1)) + ...
            Gmax_KAHP * mvec_KAHP(1,i) * (E_K-Vm_D(1,i-1)) - ...
            G_link * (Vm_D(1,i-1) - Vm_S(1,i-1)) + Iapp_D)*dt/C_D;
        Vm_D(1,i) = Vm_D(1,i-1)+dV_D;
    
        %% detect spike and increment spk_detect if detected a spike
        if spk_detect   % when Vm have dropped below down_th
            if Vm_S(1,i)>up_th 
                spk_detect=0;
                spk_count=spk_count+1;

                % When inter-spike interval is greater than the 
                % min_interval, it count as a new burst
                if (i-last_spike)>(min_interval/dt)
                    burst_count=burst_count+1;
                end
                
                last_spike = i; % update last spike
            end
        else    % when Vm have rise above up_th, not yet drop below down_th
            if Vm_S(1,i)<down_th
                spk_detect=1;
            end
        end
           
    end
    
    %% Within each trial, plot graph for the first trial, and record the 
    % number of spikes for the rest of trials to explore (G/I)-F relation
    % applied when question_number = 3,5,6
    if trial ==1 && question_number>2
        disp("When G_{link} = " +num2str(G_link*1e9) + ...
            "nS, Iapp_{dendrite} = " + num2str(Iapp_D*1e12) + "pA,"+ ...
            "and Iapp_{soma} = " + num2str(Iapp_S*1e12) + "pA")
        disp("number of spikes is "+ num2str(spk_count))
        disp("number of bursts is "+ num2str(burst_count))
        disp("average number of spikes per burst is " + num2str(spk_count/burst_count))
        figure(2)
        subplot(2,1,1)
        plot(tvec,Vm_D*1000);
        ylabel("Dendritic Vm (mV)")
        axis([0 tmax -100, 80])
        
        subplot(2,1,2)
        plot(tvec,Vm_S*1000);
        ylabel("Somatic Vm (mV)")
        xlabel("time (s)")
        axis([0 tmax -100, 80])
        
        sgtitle({"Membrane Potential for Neural Components,", " G_{link} = "+ ...
            num2str(G_link*1e9)+ "nS, Iapp_{Dendrite} = "+ ...
            num2str(Iapp_D*1e12)+ "pA, Iapp_{Soma} = "+ ...
            num2str(Iapp_S*1e12)+ "pA"}, "FontWeight", "Bold")
    end

    %% applied only for multiple trials when exploring relationships 
    % between parameters and spike number
    if trial>1        
        spk(trial-1) = spk_count;   % record number of spikes
        burst(trial-1) = burst_count;

        count=spk_count/burst_count;
        if count<5
            rate(trial-1) = count;
        else
            rate(trial-1)=1;
        end
    end
end

%% plot the relationship between G_link or Iapp_dentrite with spike number
if question_number==5
    figure(3)
    scatter(G_link_values,spk)
    title(["Correlation between G_{link} and number of spikes", ...
        "in 2-scecond simulation"])
    xlabel("G_{link} (nS)") 
    ylabel("spiking count")

    figure(4)
    scatter(G_link_values, rate)
    title(["Correlation between G_{link} and number of spikes", ...
        "per burst"])
    xlabel("G_{link} (nS)")  
    ylabel("spiking per burst")

elseif question_number==6
    figure(3)
    scatter(Iapp_values, spk)
    if mode=='s'
        title(["Correlation between somatic applied current", ...
            "and number of spikes in 2-scecond simulation"])
        xlabel("Iapp_{somatic} (pA)") 
    elseif mode=='d'
        title(["Correlation between dendritic applied current", ...
            "and number of spikes in 2-scecond simulation"])
        xlabel("Iapp_{dendritic} (pA)") 
    end
    ylabel("spiking count")

    figure(4)
    scatter(Iapp_values, rate)
    if mode=='s'
        title(["Correlation between somatic applied current", ...
        "and number of spikes per burst"])
        xlabel("Iapp_{somatic} (pA)") 
    elseif mode=='d'
        title(["Correlation between dendritic applied current", ...
        "and number of spikes per burst"])
        xlabel("Iapp_{dendritic} (pA)")  
    end
    
    ylabel("spiking per burst")
end
