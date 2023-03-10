function [fano, variance, mean] = fano_factor(spikes, dt, bin)
% Input parameters are spike vector (spikes) contains only 0 and 1, dt for 
% spike vector, and desired bin size (bin), can be a list of bin size
% values.
    fano = zeros(size(bin));
    mean = zeros(size(bin));
    variance = zeros(size(bin));

    for j = 1: length(bin)
        bin_size = bin(j);

        Nt = floor(bin_size/dt);    % number of time steps within each time bin

        Nbin = floor(length(spikes)/Nt);   % number of bins
        
        spk_num = zeros(1, Nbin);   % vector that counts number of spikes within each time bin

        for i = 1:Nbin
            spk_num(i) = length(find(spikes((i-1)*Nt+1 : i*Nt)));
        end
    
        mean(j) = sum(spk_num)/Nbin;    % mean spike time within all bins
        variance(j) = var(spk_num);     % variace of spike time
        fano(j) = variance(j)/mean(j);        % fano vector
    end

    
end