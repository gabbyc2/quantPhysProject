% Define parameters
A1 = 0.5;
A2 = 0.5;
Mu_max = 0.00024;
Mu_half = 0.00012;
M0 = 0;
Kf = 97000000; % egf association rate
Kr = 0.24;      % egf dissociation rate
Ku = 10;     % EGF/HER coupling rate
L0 = 0.000000001612; % initial EGF ligand concentration
Pr = 200000;   % cell surface egf receptors per unit cell
Ph = 60000;   % total cell surface her2 receptors per unit cell
Kc = 0.001;     % a constant

% Define initial conditions
C0 = 0; % Initial value of C = changed from 1300
T0 = 0; % Initial value of T = changed from 1300
N0 = 50000;    % Initial value of N

% Define the derivative function mu(C, T)
mu = @(C, T) (A1 * ((Mu_max * C) ./ (Mu_half + C)) + A2 * ((Mu_max * T) ./ (Mu_half + T)) + M0) * N0;

% Define the ODEs for the dynamics of C, T, and N
dCdt = @(t, C, T, N) Kf * (Pr * N - C - T) * L0 - Kr * C - Kc * (Ph * N - T) * C + Ku * T;
dTdt = @(t, C, T, N) Kc * (Ph * N - T) * C - Ku * T;
dNdt = @(t, C, T, N) mu(C, T); % Derivative of N(t) is mu(C, T)

% Define combined derivative function
odeFunc = @(t, y) [dCdt(t, y(1), y(2), y(3)); dTdt(t, y(1), y(2), y(3)); dNdt(t, y(1), y(2), y(3))];

% Solve the ODEs
[t, y] = ode45(odeFunc, [0, 100], [C0; T0; N0]);

% Extract the solution
C = y(:, 1);
T = y(:, 2);
N = y(:, 3);

% Plot N(t) over time
plot(t, N);
xlabel('Time (t)');
ylabel('Population (N)');
title('Plot of N(t) over time');