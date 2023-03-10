% Tutorial 5.1. Synaptic responses to changes in inputs
% Include main part, as well as option A and B
% Alter option at line 8 to either 0 or 1 to view the optional question
% Alter option_num at line 9 to either 'A' or 'B' to view plots
% Written by Clark Xu, last modified on 3/9/2023

clear
option=1;
option_num = 'A';

dt = 1e-4;  % 0.1 ms
tmax = 4;   % 4s simulation
tvec = 0:dt:tmax;
Nt = length(tvec);
Ntrial = 1;

if option && option_num=='A'
    Ntrial = 50;
elseif option && option_num=='B'
    Ntrial = 3;
end

% create vectors and initialize them
pre_fr_val = [20, 100, 10, 50];
pre_fr = zeros(Ntrial, Nt); % presynaptic firing rate vector
F = zeros(size(pre_fr));  % fascilitation vector
D = zeros(Ntrial, 2, Nt);  % depression vector
G = zeros(Ntrial, 3, Nt);

D(:,:,1) = 1;
F(:,1) = 1;



% generate firing rate vector
for i=1:4
    pre_fr(:,(i-1)/dt+1:i/dt) = pre_fr_val(i);
end

%  random poisson spike generated for main part and option A
random_generator = rand(size(pre_fr));
spikes = (random_generator<=dt*pre_fr); % binary list stores firing time


% parameters
delG = 1e-9;     % G increment
p0_D = 0.5;      % initial release probablity for first synaptic depression vector D1
p0_G = 0.5;     % initial release probablity for G2 and G3
p0_F = 0.2;     % initial release probablity for F
p0_D2 = 0.2;    % initial release probablity for second depression vector D2


tau_G = 0.100;      % time constant for third synaptic conductance vector
tau_D = 0.25;       % time constant for synaptic depression vector
tau_F = 0.25;       % time constant for synaptic facilitation vector

del_Gmax = 2e-9;    % del_Gmax for the second synaptic conductance vector
del_G3max = 4e-9;   % del_Gmax for the third synaptic conductance vector
f_fac = 0.25;
F_max = 1/p0_F;

if option && option_num=='A'
    Ntrial=50;
    tau_G = 0.002;  % alter time constant for G
end

% set up parameters for postsynaptic cell as a leaky integrate-and-fire neuron
% For each trial for option B, the ith trial refers to ith conductance vector
if option && option_num=='B'
    E_L = -0.065;
    E_syn = 0e-3;   % reversal potential of an excitatory synpase
    G_L = 2e-9;
    C_m = 20e-12;
    V = zeros(3,Nt);
    Vth = -0.050;
    V_reset = -0.080;    
end


%% simulaiton
for trial = 1:Ntrial        
   
    for i=2:Nt
        %% update decay of rise
        F(trial,i) = F(trial,i-1) + (1-F(trial,i-1))/tau_F*dt;    % facilitation increment    
        D(trial,:,i) = D(trial,:,i-1) + (1-D(trial,:,i-1))/tau_D*dt;  % increment Decay        
        G(trial,:,i) = G(trial,:,i-1) - G(trial,:,i-1)/tau_G*dt;    % G decays
        
        % update V for option B
        if option && option_num=='B'              
            V(trial,i) = V(trial,i-1) + (G_L*(E_L-V(trial,i-1))+ ...
                G(trial,trial,i-1)*(E_syn-V(trial,i-1)))*dt/C_m;  
            if V(trial,i)>Vth
                V(trial,i) = V_reset;
                spikes(trial,i) = 1;
            end
        end
         
        if spikes(trial,i)   % changes when fires, spikes(i) == 1, otherwise 0            
            F(trial,i) = F(trial,i) + f_fac*(F_max-F(trial,i-1));
            D(trial,1,i) = D(trial,1,i) - p0_D*D(trial,1,i);    % decrease by amount of p0*D
            D(trial,2,i) = D(trial,2,i) - p0_D2*F(trial,i-1)*D(trial,2,i-1);
            G(trial,1,i) = G(trial,1,i) + delG;     % increment by fixed delG = 1nS
            G(trial,2,i) = G(trial,2,i) + del_Gmax*p0_G*D(trial,1,i-1);   % relate to D only     
            G(trial,3,i) = G(trial,3,i) + del_G3max*p0_G*F(trial,i-1)*D(trial,2,i-1); % relate to D and F
        end
    end  
end


%% plot the graph   
if option==0 
    figure(1)       
    subplot(4,1,1)
    plot(tvec, spikes)
    for i = 1:3            
        subplot(4,1,i+1)
        Gplot = G(1,i,:);
        plot(tvec,Gplot(:));
    end
    
    figure(2)
    subplot(3,1,1)
    plot(tvec, spikes)
    for i = 1:2            
        subplot(3,1,i+1)
        Dplot = D(1,i,:);
        plot(tvec,Dplot(:));
    end
    
    figure(3)
    subplot(2,1,1)
    plot(tvec, spikes)
    subplot(2,1,2)
    plot(tvec, F)
else    % plot graphs for optional questions
    if option_num=='B'

        figure(1)
        spike_count = zeros(3,4);
        for i=1:3
            subplot(3,1,i)
            plot(tvec,V(i,:))
            
            % count the number of spikes in each second interval
            for j=1:4
                spike_count(i,j) = length(find(spikes(i,(j-1)/dt+1:j/dt)));
            end
        end
        disp(spike_count)
    elseif option_num=='A'
        Gsum = reshape(sum(G,1), 3, Nt);    % sum up 50 trials

        figure(1)
        subplot(2,1,1)
        plot(tvec,reshape(G(1,1,:),1,Nt))
        subplot(2,1,2)
        plot(tvec, Gsum(1,:))
        

        figure(2)
        subplot(2,1,1)
        plot(tvec,reshape(G(1,2,:),1,Nt))
        subplot(2,1,2)
        plot(tvec, Gsum(2,:))
    end
end


