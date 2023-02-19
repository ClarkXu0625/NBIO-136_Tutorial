% vector for different values of membrane potential
Vm_values = -0.085:0.001:0.050;

% vector of calcium concentration
Ca_conc = 0: 5e-5: 2e-3;

[alpha_m, beta_m, alpha_h, beta_h , alpha_n, beta_n] = PR_soma_gating(Vm_values);
[alpha_mca, beta_mca, alpha_kca, beta_kca, alpha_kahp, beta_kahp ] = PR_dend_gating( Vm_values, Ca_conc );

%% Set default styles for the plot
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');

%% Question 2, 
% plot all 12 of the rate constants for gating variables 
% figure 1 are the gating variables dependent on membrane potential
% figure 2 are the gating variables dependent on Calcium concentration
figure(10)
subplot(2,3,1)
plot(Vm_values,alpha_m)
hold on
plot(Vm_values,beta_m)
legend(["alpha_m", "beta_m"])

subplot(2,3,2)
plot(Vm_values,alpha_h)
hold on
plot(Vm_values,beta_h)
legend(["alpha_h", "beta_h"])

subplot(2,3,3)
plot(Vm_values,alpha_n)
hold on
plot(Vm_values,beta_n)
legend(["alpha_n", "beta_n"])

subplot(2,3,4)
plot(Vm_values, alpha_mca)
hold on
plot(Vm_values, beta_mca)
legend(["alpha_mca", "beta_mca"])

subplot(2,3,5)
plot(Vm_values, alpha_kca)
hold on
plot(Vm_values, beta_kca)
legend(["alpha_kca", "beta_kca"])


% add label 
for i=1:5
    subplot(2,3,i)
    xlabel("Membrane potential")
    ylabel("Gating variable")
end
figure(11)
plot(Ca_conc, alpha_kahp)
hold on
plot(Ca_conc, beta_kahp)
legend(["alpha_kahp", "beta_kahp"])