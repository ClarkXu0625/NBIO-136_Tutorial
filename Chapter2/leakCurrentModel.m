% Tutorial 2.1, Leaky integrate and fire model
% Written by Clark Xu

%% parameter
E_l = -0.070;
R_m = 5e6;
C_m = 2e-9;
V_reset = -0.065;
V_threshold = -0.050;
t_max = 2;
dt = 1e-4;
tau_m = C_m*R_m;
question_number=2;

% threshold for current
I_threshold = (1/R_m)*(V_threshold-E_l); 

%% Create vectors
% Io is the list of Iapp values
Io = 1:0.01:1.5;
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

%% simulation
for trial=1:Ntrial    
    firing_count = 0;   % initialize firing count to 0 every trial
    I_app = Io(trial)*ones(1,Nt);   % applied current vector
    
    %% add noise into the circuit
    switch question_number
        case 1
            noise_vec = zeros(1,Nt);    % No noise in question1
        case 2
            sigma_I =2.3e-7;
            noise_vec = rand(1,Nt)*sigma_I*sqrt(dt);
    end

    for i=2:Nt
        dV = ((E_l-V(trial,i-1))/R_m+I_app(1,i)+noise_vec(1,i))*(dt/C_m);
        V(trial,i) = V(trial,i-1)+dV;

        % Reset if exceed threshold
        if V(trial,i)>V_threshold
            V(trial,i)=V_reset;
            firing_count = firing_count+1;
        end
    
    end

    firing_rate = firing_count/2;  % simulation is 2-second long
    % firing rate limited below within 0-100Hz range
    firing_rate=min(100,firing_rate);
    firing_rate=max(0,firing_rate);
    firing(1, trial)=firing_rate;
    
    % Inter-spike interval
    % max statement is used to ensure no negative number in natural log
    ISI = tau_m*log(max(0, Io(1,trial)*R_m+E_l-V_reset))- ...
            tau_m*log(max(Io(1,trial)*R_m+E_l-V_threshold));
    % calculated firing rate
    calculated_firing(1, trial) = 1/ISI;

end

plot(Io, firing);



