function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
% Inputs are respectively: applied current vector; vector of spike times, 
% with 1 at the time of every spike; the time-step used in these vectors; 
% the time of the stimulus before a spike to begin recording; and the time 
% of the stimulus after a spike to stop recording. 
% Spike-triggered average (STA) within the given time window is returned.
    
    % set up default values for tminus and tplus if not provided.
    if (~exist('tminus', 'var'))
        tminus = 0.075;
    end

    if (~exist('tplus', 'var'))
        tplus = 0.025;
    end

    % number of time bins before and after a spike
    nminus = tminus/dt;
    nplus = tplus/dt;    
    
    % time window
    tcorr = -nminus*dt :dt: nplus*dt;

    sta = zeros(size(tcorr));   % initialize vector sta
    pks = find(spikes);         % the vector of index that simulation fires
    Nspike = length(pks);       % number of spikes
    Nignored = 0;               % number of spikes being ignored
    
    for i=1:Nspike
        if pks(i)<tminus/dt || length(Iapp)-pks(i)<tplus/dt
            Nignored = Nignored+1;
        else
            sta = sta + Iapp(pks(i)-tminus/dt : pks(i)+tplus/dt);
        end
    end
    
    % take the average
    sta = sta/(Nspike-Nignored);
end