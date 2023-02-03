dt = 0.001;
t_max=10;
t = 0:dt:t_max;
Nt=length(t);
Nparameter = 2;

x = zeros(1, Nt);
v = zeros(1, Nt);



%% parameters
w = [1 2];
x0 = 1;
v0 = 1;

x(1, 1) = x0;
v(1, 1) = v0;

for parameter = 1:Nparameter
    for i =2:Nt
        dv = (-w(parameter)^2*x(1,i-1))*dt;
        dx = v(1,i-1)*dt;
        v(1, i) = v(1, i-1)+dv;
        x(1, i) = x(1, i-1)+dx;
        
    end
    plot(t, x);
    hold on
    
end

xlabel("time")
ylabel("distance")