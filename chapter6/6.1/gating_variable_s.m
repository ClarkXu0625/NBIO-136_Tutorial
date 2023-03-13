function s_curve = gating_variable_s(r_value, p_r, alpha_0, tau_s, tau_D)

    d_curve = ones(1,length(r_value));
    s_curve = zeros(1,length(r_value));

    if (exist('tau_D'))
        d_curve = 1. / (p_r*r_value*tau_D + 1);

    end
    

    for i = 1:length(r_value)
        s_curve(i) = alpha_0*d_curve(i)*r_value(i)*tau_s / ...
            (1 + alpha_0*d_curve(i)*r_value(i)*tau_s);
    end