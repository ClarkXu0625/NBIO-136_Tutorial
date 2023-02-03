% This is code written for tutorial 2.2, Modeling the refractory period
% Written by Clark Xu

%% Parameter
E_l = -0.070;
R_m = 1e8;
C_m = 1e-10;
t_max = 2;
dt = 1e-4;
t = 0:dt:t_max;
Nt = length(t);

question_number = 1;

%% setting up Vectors
Vm = zeros(1, Nt);




switch question_number
    case 1
        V_th = -0.050;
        V_reset = -0.065;
        tau_ref = 2.5e-3;   % tau refractory
    
    
end    


Iapp = 600e-12;
for i = 2:Nt


    switch question_number
        case 1
            dV = ((E_l-Vm(1,i-1))/R_m+Iapp)*(dt/C_m);
            Vm(1,i)=Vm(1,i-1)+dV;
            if Vm(1,i)>V_th
                Vm(1,i)=V_th;
            end
            
    end

    

end

plot(t, Vm)