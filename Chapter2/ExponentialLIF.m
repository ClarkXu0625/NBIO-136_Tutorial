% Tutorial 2.3, Models based on extensions of the LIF neuron.
% Question 1 simulates an LIF model, and question 2 simulates AELIF model
% Written By Clark Xu, Feb. 7, 2023
%%%%%%%%%%%%

% Accepted inputs are 1, 2, 3, and 4. 1 refers to question 1(a), 2 refers 
% to question 1(b), 3 refers to 2(a), 4 refers to 2(b)
question_number=3;

%% Parameters
E_l = -0.075;   % reverse potential
Vth = -0.050;   % threashold for membrane potential
Vreset = -0.080;    % membrane potential for reset
switch question_number
    case {1,2}        
        R_m = 1e8;
        C_m = 100e-12;
        E_K = -0.080;
        delta_G_RSA = 1e-9;
        tau_RSA = 0.200;
        tau_m = C_m*R_m;
    case {3,4}
        Vmax = 0.100;
        delta_th=0.002;
        G_L = 10e-9;
        C_m = 100e-12;
        a = 2e-9;
        b = 0.02e-9;
        tau_SRA = 0.200;
end


%% setting up vectors
switch question_number
    case {1,3}
        tmax = 1.5;
    case {2,4}
        tmax = 5;
end

dt = 0.0001;    % time step is 0.1 ms
t = 0:dt:tmax;  % t vector
Nt = length(t);
Iapp_values = (100:15:550)*1e-12;   % vector of Iapp values used in 1(b) and 2(b)

% Set Iapp.
% Ntrials is the number of different applied currents simulation
switch question_number
    case {1,3}  % Quetion 1 and 3 simulates step pulse with only one Iapp value
        Ntrial=1;
        Iapp_vector=zeros(1,Nt);
        Iapp_vector(0.5/dt:1/dt) = 500e-12; % step pulse
    case {2,4}  % Qeustion 2 and 4 simulates multiple contant Iapp value
        Ntrial=length(Iapp_values);
        Iapp_vector=ones(Ntrial,Nt);        
end

% for each trial, count the ISI for each spike by recording the
% time which spikes occur; first space records the initial INI, and the
% second space records the steady state INI.
ISI = zeros(Ntrial, 2);

% Set up vectors for Vm, G_SRA, and I_SRA
V = zeros(Ntrial, Nt);  % Membrane Potential
V(:,1)=E_l; % initialize membrane potential
G_SRA = zeros(Ntrial,Nt);   % conductance SRA
I_SRA = zeros(Ntrial,Nt);

%% simulation
for trial=1:Ntrial
    
    % set Iapp to constant values
    switch question_number
        case {2,4}
            Iapp_vector(trial,:) = Iapp_values(1,trial)*Iapp_vector(trial,:);
    end
    

    last_spike = 0; %initialize time for last spike every trial
    interval=0; % initialize INI every trial
    initial_interval = 0;

    for i = 2:Nt
        
        switch question_number
            case {1,2}
                dG_SRA = -G_SRA(trial,i-1)*dt/tau_RSA;
                G_SRA(trial,i) = G_SRA(trial,i-1)+dG_SRA;
            
                dV = ((E_l-V(trial,i-1))/R_m+G_SRA(trial,i)*(E_K-V(trial,i-1))+ ...
                    Iapp_vector(trial,i))*dt/C_m;
                V(trial,i)=V(trial,i-1)+dV;
            case {3,4}
                dI_SRA = (a*(V(trial,i-1)-E_l)-I_SRA(trial,i-1))*dt/tau_SRA;
                I_SRA(trial,i) = I_SRA(trial,i-1)+dI_SRA;

                dV = (G_L*(E_l-V(trial,i-1)+delta_th*exp((V(trial,i-1)-Vth)/delta_th))- ...
                    I_SRA(trial,i)+Iapp_vector(trial,i))*dt/C_m;
                V(trial,i) = V(trial,i-1)+dV;
                
        end

        % spike occurs
        if V(trial,i)>Vth
            V(trial,i)=Vreset;  % reset membrane potential

            switch question_number
                case {1,2}  % increment G_SRA
                    G_SRA(trial,i)=G_SRA(trial,i)+delta_G_RSA;
                case {3,4}  % increment I_SRA
                    I_SRA(trial,i)=I_SRA(trial,i)+b;
            end
            
            interval=i*dt-last_spike;   % INI updated every spike

            if last_spike==0
                ISI(trial,1)=interval;  % fist inter-spike interval      
                initial_interval = interval;
            end  

            last_spike=i*dt;    % update last spike           
        end
        
        ISI(trial,2)=interval;  % last interval is the steady-state ISI

        % if last spike is not lasting until the end of the simulation, the
        % system did not enter a steady state. The steady state interval
        % will be 0
        if last_spike<4
            ISI(trial,2) = 0;
        end
        
    end
end

% Get the inverse of inter-spike interval value as the firing rate.
% replace all Inf from zero division in matrix by 0
firing_rate=1./ISI;
firing_rate(isinf(firing_rate)) = 0;    

%% plot graphs
clf     % clear the figure

% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

switch question_number
    case 1
        figure(1)
        subplot(3,1,1);
        plot(t, Iapp_vector);
        ylabel("Iapp (A)")
        subplot(3,1,2);
        plot(t, V);
        ylabel("membrane potential (V)")
        subplot(3,1,3);

        plot(t, G_SRA);
        ylabel("adaptation conductance (F)")
        xlabel("time (s)")
        
    case 3
        figure(1)
        subplot(2,1,1);
        plot(t, Iapp_vector);
        ylabel("Iapp (A)")
        subplot(2,1,2);
        plot(t, V);
        ylabel("membrane potential (V)")
        xlabel("time(s)")
        
    case {2,4}
        figure(1)
        ylabel("inverse of ISI")
        scatter(Iapp_values,firing_rate(:,1),'filled','<');
        hold on
        plot(Iapp_values,firing_rate(:,2));
        ylabel("firing rate (Hz)")
        legend(["initial state", "steady state"], 'Location', 'northwest');
        xlabel("Iapp (A)")
        ylabel("Spike Rate (Hz)")
        yline(50, "--", "firing rate = 50Hz" ,'HandleVisibility','off')
end

