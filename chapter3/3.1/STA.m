function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
    
    if (~exist('tminus', 'var'))
        tminus = 0.075;
        nminus = tminus/dt;
    else
        nminus = tplus/dt;
    end

    if (~exist('tplus', 'var'))
        tplus = 0.025;
        nplus = tplus/dt;
    else
        nplus = tplus/dt;
    end

    % number of time bins before and after a spike
    %nminus = tminus/dt;
    %nplus = tplus/dt;    
    
    % time window
    tcorr = -nminus*dt :dt: nplus*dt;

    sta = zeros(size(tcorr));
    pks = find(spikes);
    Nspike = length(pks);
    Nignored = 0;
    
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