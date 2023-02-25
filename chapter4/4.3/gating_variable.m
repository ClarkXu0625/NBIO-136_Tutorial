%% Output is the change in gating variable. Take alpha, beta, gating 
%  variable at the last time step, and time step size input.
function d = gating_variable(alpha, beta, previous, dt)
    eq = alpha/(alpha+beta);    % steady state
    tao = 1/(alpha+beta);   % time constant
    d = (eq-previous)*dt/tao;
end