clc; clear;

% Constants from the table
kf = 9.7e7; % EGF association rate constant (M^−1 min^−1)
kr = 0.24; % Rate constant for EGF-EGFR and HER2 coupling (min^−1)
ku = 1e-3; % Rate constant for EGF-EGFR-HER2 dissociation ((#cell)^−1 min^−1)
kc = 0.1; % can be adjusted
kL = 10; % Rate constant for EGF degradation (min^−1)
L0 = 1.612e-9; % Initial EGF ligand concentration (M)
N0 = 0.5e6; % Initial cell density (#/volume)
muL0 = 0.00; % Cell proliferation rate constant independent of EGFR and HER2 (min^-1)
lambdaL = 0.0; % Ligand degradation rate (Molecule/min)
A1 = 0.5; % nonnegative constant
A2 = 0.5; % nonnegative constant
mu1max = 0.0143; % maximum specific cell growth rate constant (time−1)
mu2max = 0.0143; % maximum specific cell growth rate constant (time−1)
mu1half = 0.00715; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu2half = 0.00715; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu0 = 0.01; % cell proliferation rate independent of EGFR and HER2 receptors (time−1)
rho_R = 0.1; % can be adjusted
rho_H = 0.1; % can be adjusted
S_L = 0.1; % can be adjusted
lambda_d = 0.1; % can be adjusted

% Time vector
t = linspace(0,100,10000); 

% Initial conditions
C = 200;
T = 300;
L = L0;


for i = 2:length(t)
    % Calculate mu_CT
    mu_CT = A1 * mu1max * C / (mu1half + C) + A2 * mu2max * T / (mu2half + T) + mu0;
end

% Cell density over time
N_t = N0*exp(mu_CT.*t);

mu_Ct = 0;

% Differential equations for C, T, and L
for i = 2:length(t)
    % Calculate mu_CT
    mu_CT = A1 * mu1max * C / (mu1half + C) + A2 * mu2max * T / (mu2half + T) + mu0;
    
    dC_dt = kf*(rho_R*N_t(i-1) - C - T)*L - kr*C - kc*(rho_H*N_t(i-1) - T)*C + ku*T;
    dT_dt = kc*(rho_H*N_t(i-1) - T)*C - ku*T;
    dL_dt = -kf*(rho_R*N_t(i-1) - C - T)*L + kr*C + S_L - lambda_d*L;
    
    % Update C, T, and L
    C = C + dC_dt;
    T = T + dT_dt;
    L = L + dL_dt;
end

% Plotting the graph for cell density over time
figure;
plot(t,N_t,'LineWidth',2);
xlabel('Time');
ylabel('Cell Density');
title('Cell Density Over Time');
grid on;
