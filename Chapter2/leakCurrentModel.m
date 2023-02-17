% Tutorial 2.1, Leaky integrate and fire model
% Written by Clark Xu
% Feb. 3, 2023
% Last modified on Feb. 14, 2023

%% parameter
E_l = -0.070;
R_m = 5e6;
C_m = 2e-9;
V_reset = -0.065;
V_threshold = -0.050; 
t_max = 2;  % 2 second simulation
dt = 1e-4;  % time scale is 0.1 ms
tau_m = C_m*R_m;
question_number=2; % Accepted inputs are 1 and 2 for question1 and 2
compare = 1;    % 1 if want to compared counted and calculated rates (1d), 0 otherwise (1c)

% threshold for current
I_threshold = (1/R_m)*(V_threshold-E_l);

%% Create vectors
% Io is the list of Iapp values
Io = 0.8:0.01:1.6;
Ntrial=length(Io);
Io = I_threshold*Io;

t = 0:dt:t_max;
Nt = length(t);
V = zeros(Ntrial, Nt);
V(:,1) = E_l;   % initialize membrane potential by E_l

% firing rate counted from plotting
firing = zeros(1, Ntrial);
% firing rate calculated from inverse of ISI
calculated_firing = zeros(1, Ntrial);

switch question_number
    case 1
        sigma_num = 1;
    case 2
        sigma = [0, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 2e-2];
        sigma_num = length(sigma);
end



% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%% simulation
for num = 1:sigma_num
    for trial=1:Ntrial    
        firing_count = 0;   % initialize firing count to 0 every trial
        I_app = Io(trial)*ones(1,Nt);   % applied current vector
        
        %% add noise into the circuit
        switch question_number
            case 1
                noise_vec = zeros(1,Nt);    % No noise in question1
            case 2
                sigma_I =sigma(num);  % choose one sigma value from list
                noise_vec = randn(1,Nt)*sigma_I*sqrt(dt);
        end
    
        for i=2:Nt
            dV = ((E_l-V(trial,i-1))/R_m+I_app(1,i))*(dt/C_m);
            V(trial,i) = V(trial,i-1)+dV+noise_vec(1,i);
    
            % Reset if exceed threshold
            if V(trial,i)>V_threshold
                V(trial,i)=V_reset;
                firing_count = firing_count+1;
            end
        
        end
    
        firing_rate = firing_count/2;  % simulation is 2-second long
        % firing rate limited below within 0-100Hz range
        firing_rate=max(0,firing_rate);
        firing(1, trial)=firing_rate;
        
        % Inter-spike interval
        % max statement is used to ensure no negative number in natural log
        ISI = tau_m*log(max(0, Io(1,trial)*R_m+E_l-V_reset))- ...
                tau_m*log(max(0, Io(1,trial)*R_m+E_l-V_threshold));
        
        % calculated firing rate
        calculated_firing(1, trial) = 1/ISI;
    
    end

    switch question_number
    case 1
        % plot membrane potential
        figure(2)
        plot(0:dt:0.2, V(22, 1:0.2/dt+1))
        xlabel("time (s)")
        ylabel("V membrane (V)s")
        title("Membrane potential when applied current is "+string(Io(22)) + "A")

        % plot I-f graph
        figure(1)
        plot(Io, firing);
        xlabel("Applied Current");
        ylabel("firing rate(Hz)");
        hold on
        if compare
            plot(Io, calculated_firing)            
            legend("counted", "calculated", "Location","northwest")            
        end
        
        case 2     % question 2   
        plot(Io, firing);
        xlabel("Applied Currentm (A)")
        ylabel("firing rate (Hz)")
        title("question 2(b)")
        hold on        
    end
end

% legend setup for plots in question 2
if question_number==2
    legend_names = "sigma_I = " + string(sigma(1:sigma_num));
    legend(legend_names, "Location","northwest")
end






