% Tutorial 2.1, Leaky integrate and fire model
% Written by Clark Xu

%% parameter
E_l = -0.070;
R_m = 5e6;
C_m = 2e-9;
V_reset = -0.065;
V_threshold = -0.050;
t_max = 2;
dt = 0.0001;
Ntrial=100;
tau_m = C_m*R_m;

%% Create vectors
t = 0:dt:t_max;
Nt = length(t);
V = zeros(Ntrial, Nt);


V(:,1) = E_l;


I_0 = (1/R_m)*(V_threshold-E_l); %threshold

% Io is the list of Iapp values
Io = rand(1, Ntrial)+1;
disp(Io)
Io = I_0*Io;
disp(Io(1,1))
firing = zeros(1, Ntrial);
calculated_firing = zeros(1, Ntrial);


for trial=1:Ntrial

    firing_count = 0;
    I_app = Io(trial)*ones(1,Nt);
    
    for i=2:Nt
        dV = ((E_l-V(trial,i-1))/R_m+I_app(1,i))*(dt/C_m);
        V(trial,i) = V(1,i-1)+dV;
        if V(trial,i)>V_threshold
            V(trial,i)=V_reset;
            firing_count = firing_count+1;
        end
    
    end

    % firing rate limited below 100 Hz
    firing_count=min(100,firing_count/2);
    firing(1, trial)=firing_count;
    
    % calculated firing rate
    calculated_firing(1, trial) = 1/(tau_m*log(Io(1,trial)*R_m+E_l-V_reset)- ...
        tau_m*log(Io(1,trial)*R_m+E_l-V_threshold));

end

scatter(calculated_firing, firing);



