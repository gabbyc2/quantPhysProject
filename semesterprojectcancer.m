% Constants
A1 = 0.5; % nonnegative constant
A2 = 0.5; % nonnegative constant
mu1max = 0.7; % maximum specific cell growth rate constant (time−1)
mu2max = 0.6; % maximum specific cell growth rate constant (time−1)
mu1half = 100; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu2half = 150; % number of occupied receptors required to generate a half-maximal response (#/unit volume)
mu0 = 0.01; % cell proliferation rate independent of EGFR and HER2 receptors (time−1)

% Time vector
t = linspace(0,10,100); 

% C and T values can be modified as per requirement 
C = 200;
T = 300;

% Calculating mu(C,T) over time using the given equation
mu_CT = A1 * mu1max * C ./ (mu1half + C) + A2 * mu2max * T ./ (mu2half + T) + mu0;

% Cell density over time assuming an initial density of cells N(0)=100 
N_0=100;
N_t=N_0*exp(mu_CT.*t);

% Plotting the graph for cell density over time
figure;
plot(t,N_t,'LineWidth',2);
xlabel('Time (hrs)');
ylabel('Cell Density (# ');
title('Cell Density Over Time');
grid on;