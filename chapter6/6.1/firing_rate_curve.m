function f = firing_rate_curve(r_0, r_max, x, sigma, S)
% input is r_0, r_max, sigma, and S.
% output firing rate curve would be the same size as input S
    N = length(S);

    f = zeros(1,N);

    for i = 1:N
        S_value = S(i);
        if S_value > 0
            f_value = r_0 + r_max*((S_value^x)/(S_value^x+sigma^x));
        else
            f_value = r_0;
        end
        f(i) = f_value;
    end

end