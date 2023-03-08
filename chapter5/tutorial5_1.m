% Tutorial 5.1. Synaptic responses to changes in inputs


dt = 1e-4;  % 0.1 ms
tmax = 4;
tvec = 0:dt:tmax;
Nt = length(tvec);

pre_fr_val = [20, 100, 10, 50];
pre_fr = zeros(size(tvec)); % presynaptic firing rate vector
D = zeros(size(tvec));  % depression vector
F = zeros(size(tvec));  % fascilitation vector
G = zeros(size(tvec));
G2 = zeros(size(tvec));

% generate firing rate vector
for i=1:4
    pre_fr((i-1)/dt+1 : i/dt) = pre_fr_val(i);
    disp(pre_fr_val(i))
end

% poisson spike
random_generator = rand(size(tvec));
spikes = (random_generator<=dt*pre_fr); % binary list stores firing time
spike_index = find(spikes);     % list of spiking time


% parameters
delG = 1e-9;    % G increment
p0_D = 0.5;       % initial release probablity
p0_G = 0.5;
p0_F = 0.2;
tau_G = 0.100;  % G_syn time constant
tau_D = 0.25;
tau_G2 = 0.100;
tau_F = 0.25;
del_Gmax = 2e-9;
f_fac = 0.25;
F_max = 1/p0_F;



%% simulaiton
for i=2:Nt
    G(i) = G(i-1) - G(i-1)/tau_G*dt;    % G decays
    D(i) = D(i-1) + (1-D(i-1))/tau_D*dt;  % increment Decay
    G2(i) = G2(i-1) - G2(i-1)/tau_G2*dt;    % G decays
    F(i) = F(i-1) + (1-F(i-1))/tau_F*dt;

    if ismember(i,spike_index)>0   % increment when fires
        G(i) = G(i) + delG;
        D(i) = D(i) - p0_D*F(i)*D(i);
        G2(i) = G2(i-1) + del_Gmax*p0_G*D(i-1);
        F(i) = F(i) + f_fac*(F_max-F(i-1));
    end


end

subplot(2,1,1)
plot(tvec,G);
subplot(2,1,2)
plot(tvec, G2);
