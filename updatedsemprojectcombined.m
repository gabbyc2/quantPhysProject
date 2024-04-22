clc; clear;

%mu1max and m1half messed with

% Constants from the table
kf = 9.7e7; % EGF association rate constant (M^1 min^1)
kr = 0.24; % Rate constant for EGF-EGFR and HER2 coupling (min^1)
ku1 = 500; % Rate constant for EGF-EGFR-HER2 dissociation ((#cell)^1 min^1)
kc1 = 500000; % can be adjusted
ku2 = 1e-7; % Rate constant for EGF-EGFR-HER2 dissociation ((#cell)^1 min^1)
kc2 = 0.9; % can be adjusted
kL = 10; % Rate constant for EGF degradation (min^1)
L0 = 1.612e-9; % Initial EGF ligand concentration (M)
N0 = 0.5e6; % Initial cell density (#/volume)
muL0 = 0.00; % Cell proliferation rate constant independent of EGFR and HER2 (min^-1)
lambdaL = 0.0; % Ligand degradation rate (Molecule/min)
A1 = 0.5; % nonnegative constant
A2 = 0.5; % nonnegative constant
mu1max = 0.0024; % maximum specific cell growth rate constant (time1)
mu2max = 0.0024; % maximum specific cell growth rate constant (time1)
mu1half = 1300; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu2half = 1300; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu0 = 0.01; % cell proliferation rate independent of EGFR and HER2 receptors (time1)
rho_R = 10000; % can be adjusted
rho_H1 = 600000; % overexpression
rho_H2 = 10000;% normal expression
S_L = 0.1; % can be adjusted
lambda_d = 0.1; % can be adjusted
N_max = 1e8;

% Time vector
t = linspace(0,100,10000); 

% Initial conditions
C1 = 4200;%over
T1 = 4300;%over
C2 = 200;%normal
T2 = 300;%normal
L1 = L0;%over
L2 = L0;%normal


for i = 2:length(t)
    % Calculate mu_CT_Overexpression
    mu_CT1 = A1 * mu1max * C1 / (mu1half + C1) + A2 * mu2max * T1 / (mu2half + T1) + mu0;
    % Calculate mu_CT_normal expression
    mu_CT2 = A1 * mu1max * C2 / (mu1half + C2) + A2 * mu2max * T2/ (mu2half + T2) + mu0;
end

% Cell density over time
N_t_over = (N0*exp(mu_CT1.*t)); %overexpression
N_t_normal = (N0*exp(mu_CT2.*t)); 

g_n = mu1max * (1 - N_t_normal(i-1) / N_max);
N_t_new = N_t_over + g_n;

mu_Ct1 = 0;
mu_CT2 = 0;

% Differential equations for C, T, and L
for i = 2:length(t)
    % Calculate mu_CT
    mu_CT1 = A1 * mu1max * C1 / (mu1half + C1) + A2 * mu2max * T1 / (mu2half + T1) + mu0;
    mu_CT2 = A1 * mu1max * C2 / (mu1half + C2) + A2 * mu2max * T2 / (mu2half + T2) + mu0;
    
    dC_dt1 = kf*(rho_R*N_t_over(i-1) - C1 - T1)*L1 - kr*C1 - kc1*(rho_H1*N_t_over(i-1) - T1)*C1 + ku1*T1;%over
    dC_dt2 = kf*(rho_R*N_t_normal(i-1) - C2 - T2)*L2 - kr*C2 - kc2*(rho_H2*N_t_normal(i-1) - T2)*C2 + ku2*T2;%normal
    dT_dt1 = kc1*(rho_H1*N_t_over(i-1) - T1)*C1 - ku1*T1; %over
    dT_dt2 = kc2*(rho_H2*N_t_normal(i-1) - T2)*C2 - ku2*T2;%normal
    dL_dt1 = -kf*(rho_R*N_t_over(i-1) - C1 - T1)*L1 + kr*C1 + S_L - lambda_d*L1;%over
    dL_dt2 = -kf*(rho_R*N_t_normal(i-1) - C2 - T2)*L2 + kr*C2 + S_L - lambda_d*L2;%normal
    
    % Update C, T, and L
    C1 = C2 + dC_dt1; %over
    T1 = T1 + dT_dt1;%over
    C2 = C2 + dC_dt1;%normal
    T2 = T2 + dT_dt2;%normal
    L1 = L1 + dL_dt1;%over
    L2 = L2 + dL_dt2;%normal
end

% Plotting the graph for cell density over time
figure;

plot(t,N_t_normal,'LineWidth',2);
hold on
plot(t,N_t_new,'LineWidth',2);
xlabel('Time (Hours)');
ylabel('Number of Cells');
title('Cell Population Over Time');
grid on;
legend('N_t_new','N_t_over')