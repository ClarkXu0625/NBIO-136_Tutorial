% This is code written for NBIO-136 couse
% Tutorial 2.2, Modeling the refractory period in a leaky current model.
% This code plots the f-I curve, apply one of following voltage-update rule
% 1) Forced Voltage Clamp
% 2) Threshold Increase
% 3) Refractory Conductance with threshold increase
% Alter variable question_number from either 1, 2, or 3 to choose one of
% the rules listed above.
% Written by Clark Xu, on Feb. 4, 2023.

%% Parameter
E_l = -0.070;
R_m = 1e8;
C_m = 1e-10;
t_max = 2;
dt = 1e-4;
t = 0:dt:t_max;
Nt = length(t);

question_number = 1;    % Allowed inputs are [1, 2, 3]

%% Set up Vectors
Vm = zeros(1, Nt);
Iapp_vector=(220:1:600)*1e-12;  % applied current varies from 220 to 600 pA
Ntrial=length(Iapp_vector);
rates=zeros(1,Ntrial);

%% Set up basic conditions for each question
switch question_number
    case 1
        V_th = -0.050;  % Membrane potential threshold
        V_reset = -0.065;
        tau_ref = 2.5e-3;   % refractory period is 2.5ms
        last_spike_time = -tau_ref;

    case {2, 3}
        % condition in question number 2
        V_reset = -0.065;
        V0_th = -0.050; % initial condition of threshold
        Vtheshold_vector=zeros(1,Nt);   % Vector of Vth
        Vtheshold_vector(1,1) = V0_th;  % initialize theshold to initial condition
        tau_Vth = 1e-3; % refractory time constant
        Vth_max = 0.200;    % maximum threshold
    
        % condition for question 3 only
        E_K = -0.080;
        G_ref_vector=zeros(1,Nt);   % vector for refractory conductance
        tau_Gref = 0.2e-3;
        dG = 2e-6;
end    


for trial=1:Ntrial
    Iapp = Iapp_vector(trial);
    firing_count=0;
    %% simulation start here
    for i = 2:Nt
        % question 2 and 3
        if question_number~=1
                dVth = (V0_th-Vtheshold_vector(1,i-1))*(dt/tau_Vth);
                Vtheshold_vector(1,i)=Vtheshold_vector(1,i-1)+dVth;
            if question_number==3
                dG_ref = -G_ref_vector(1,i-1)*(dt/tau_Gref);
                G_ref_vector(1,i)=G_ref_vector(1,i-1)+dG_ref;
            end          
        end
    
        
        dV = ((E_l-Vm(1,i-1))/R_m+Iapp)*(dt/C_m);
    
        % add refractory conductace element to dV
        if question_number==3
            dV = dV + G_ref_vector(1,i)*(E_K-Vm(1,i-1))*(dt/C_m);
        end
    
        Vm(1,i)=Vm(1,i-1)+dV;
    
        if question_number==1       
            % Reset the membrane potential when exceeding the threshold;
            % set the time of last spike to now.
            % Increment firing count
            if Vm(1,i)>V_th
                firing_count = firing_count+1;
                Vm(1,i) = V_reset;
                last_spike_time = i*dt;
            end
    
            % set membrane potential to Vreset for a duration of refractory
            % period after a spike.
            if (last_spike_time+tau_ref)>(i*dt)
                Vm(1,i)=V_reset;
            end
            
    
        else
            % Reset the membrane potential when exceeding the threshold;
            % change threshold of membrane potential to maximum level
            % Increment firing_count
            if Vm(1,i)>Vtheshold_vector(1,i)
                
                firing_count=firing_count+1;
                Vtheshold_vector(1,i)=Vth_max;
    
                % increase refractory conductance after a spike for question 3
                if question_number==3
                    G_ref_vector(1,i)=G_ref_vector(1,i)+dG;
                else
                    % Reset membrane potential in question 2
                    % Not resetting the membrane potential for question 3    
                    Vm(1,i)=V_reset;
                end
            end        
        end
    
    end  
    firing_rate=firing_count/t_max;
    rates(1,trial)=firing_rate;

end

plot(Iapp_vector,rates);


