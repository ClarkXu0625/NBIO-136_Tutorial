% Tutorial 4.1, the Hodgkin-Huxley model as an oscillator.
% For each question number, it matches one of questions a-f in the
% tutorial. Each question has different applied current and initial
% condition, but all question share the parameters and model.
% By Clark Xu, Feb. 10, 2023

% accepted inputs for variable question_number are [1, 2, 3, 4, 5, 6]
% Each input matches question number [a, b, c, d, e, f]
question_number=3; 


%% parameter
G_leak = 30e-9; % leak conductance
Gmax_Na = 12e-6;    % maximum sodium conductance
Gmax_K = 3.6e-6;    % maximum delayed rectifier condunctance
E_Na = 0.045;   % sodium reversal potential
E_K = -0.082;   % potassium reversal potential
E_l = -0.060;   % leak reversal potential
C_m = 100e-12;  % membrane capacitance


%% time vector
tmax = 0.35;    % 350 ms simulation
dt = 0.02e-6;  % 0.02us
tvec = 0:dt:tmax;
Nt = length(tvec);
Ntrial=1;

%% vector for Iapp
start_time=0;   % start time of the first applied current pulse
Vm0 = 0;    % initial membrane potential

% each of following gating variable initate with 0, unless specified
m0 = 0; 
n0 = 0;
h0 = 0;

% switch different applied current for each question
switch question_number
    case 1
        % Set up parameters for each question.  
        % For question(a), no applied current                 
        duration = 0;   % duration of applied current
        repetition = 0; % number of pulses
        amplitude = 0;  % amplitude of applied current pulses
        Ibaseline = 0;  % baseline of applied current
        delay = 0;  % delay time
    case 2
        % question(b)        
        start_time=0.1;
        duration = 0.1;
        repetition = 1;
        amplitude = 0.22e-9;
        Ibaseline = 0;
        delay = 0;
    case 3
        % question(c)
        duration = 5e-3;
        delay = [5e-3, 10e-3, 15e-3, 20e-3, 25e-3];
        repetition = 10;
        amplitude = 0.22e-9;
        Ibaseline = 0;
    case 4
        % question(d)
        duration = 5e-3;    % 5ms pulse
        delay = 20e-3;  % 20ms delay
        repetition = 10;    % 10 repetitions
        amplitude = 0; % inhibitory pulse is 0
        Ibaseline = 0.6e-9;
        Vm0 = -0.065;
        m0 = 0.05;
        n0 = 0.35;
        h0 = 0.5;
    case 5
        % question(e)
        start_time = 0.1;
        duration = 5e-3;
        delay = 0;
        repetition = 1;
        amplitude = 1e-9;
        Ibaseline = 0.65e-9;
        Vm0 = -0.065;
        m0 = 0.05;
        n0 = 0.35;
        h0 = 0.5;
    case 6
        % question(f)
        start_time=0.1;
        duration = 5e-3;
        delay = 0;
        repetition = 1;
        amplitude = 1e-9;
        Ibaseline = 0.7e-9;
        Vm0 = -0.065;
end

%% vector for membrane potentials, m, h, and n
Vm = zeros(Ntrial, Nt);
mvec = zeros(size(Vm));
hvec = zeros(size(Vm));
nvec = zeros(size(Vm));
Iapp = ones(size(Vm));

%% set up Iapp and initial conditions for Vm
i=1;
j=start_time/dt+1;
Iapp = Iapp*Ibaseline;
while i<=repetition
    Iapp(1,j:j+duration/dt) = amplitude;
    i = i+1;
    j = j + (duration+delay)/dt;
end  
Vm(:,1) = Vm0;

%% set initial conditions for gating variable if specified in the question
mvec(:,1) = m0;
nvec(:,1) = n0;
hvec(:,1) = h0;


%% simulation starts here
for i = 2:Nt
    % gating variable m
    alphaM = (100000 * (-Vm(1,i-1)-0.045))/ ...
        (exp(100 * (-Vm(1,i-1)-0.045))-1);
    betaM = 4000 * exp((-Vm(1,i-1)-0.070)/0.018);
    dm = (alphaM*(1-mvec(1,i-1)) - betaM*mvec(1,i-1))*dt;
    mvec(1,i) = mvec(1,i-1)+dm;

    % gating variable h
    alphaH = 70 * exp(50*(-Vm(1,i-1)-0.070));    
    betaH = 1000/(1+exp(100*(-Vm(1,i-1)-0.040)));
    dh = (alphaH*(1-hvec(1,i-1)) - betaH*hvec(1,i-1))*dt;
    hvec(1,i) = hvec(1,i-1)+dh;

    % gating variable n
    alphaN = 10000 * (-Vm(1,i-1)-0.060)/ ...
        (exp(100*(-Vm(1,i-1)-0.060))-1);
    betaN = 125 * exp((-Vm(1,i-1)-0.070)/0.08);
    dn = (alphaN*(1-nvec(1,i-1)) - betaN*nvec(1,i-1))*dt;
    nvec(1,i) = nvec(1,i-1) + dn;

    % membrane potential
    dVm = G_leak * (E_l-Vm(1,i-1)) + ...
        Gmax_Na * (mvec(1,i-1)^3) * hvec(1,i-1) * (E_Na-Vm(1,i-1)) + ...
        Gmax_K * (nvec(1,i-1)^4) * (E_K - Vm(1,i-1)) + Iapp(1,i-1);
    Vm(1,i) = Vm(1,i-1) + dVm*(dt/C_m);
end

figure(1)
subplot(2,1,1);
plot(tvec,Iapp);
ylabel("applied current")
subplot(2,1,2);
plot(tvec, Vm);
ylabel("membrane potential")
xlabel("time")

figure(2)
subplot(3,1,1)
plot(tvec,mvec)
ylabel("gating variable m")
subplot(3,1,2)
plot(tvec,nvec)
ylabel("gating variable n")
subplot(3,1,3)
plot(tvec,hvec)
ylabel("gating variable h")
xlabel("time")