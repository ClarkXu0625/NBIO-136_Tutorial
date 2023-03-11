function f = firing_rate_curve(r_0, r_max, sigma, S)
% input is r_0, r_max, sigma, and S.
% output firing rate curve would be the same size as input S
    N = length(S);
    f = zeros(1,N);
    for i = 1:N
        if S>0
            f(i) = r_0 + r_max*(S^x/(S^x+sigma^x));
        else
            f(i) = r_0;
        end
    end

end