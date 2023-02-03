%% parameter
alpha = [1 2 4];
beta = [2 1 1];
Nparameter=3;

S0 = [0.25 0.5 0.75];


dt = 0.001;
t_max=10;
t = 0:dt:t_max;
Nt=length(t);

S = zeros(1, Nt);
S(1,1)=S0(2);

for parameter= 1:Nparameter
    S = zeros(1, Nt);
    S(1,1)=S0(parameter);
    for i = 2: Nt
        dS = (alpha(1)*(1-S(1,i-1))-beta(1)*S(1,i-1))*dt;
        S(1,i)=S(1,i-1)+dS;
       
    end
    plot(t, S);
    hold on
end


