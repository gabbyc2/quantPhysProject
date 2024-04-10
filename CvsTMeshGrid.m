%parameters
A1 = .5;
A2 = .5;
M0 = 0;
Mu_max = 0.00024; %cells per min
Mu_half = 0.00012; %cells/min 
% Define function mu(C, T)
mu = @(C, T) A1 * ((Mu_max * C) ./ (Mu_half + C)) + A2 * ((Mu_max * T) ./ (Mu_half + T)) + M0;

% Define range for C and T
C_range = linspace(1300, 4800, 100);
T_range = linspace(1300, 4800, 100);

% Create grid of C and T values
[C, T] = meshgrid(C_range, T_range);

% Evaluate mu(C, T) for all combinations of C and T
mu_values = mu(C, T);

% Plot the result
surf(C, T, mu_values);
xlabel('C');
ylabel('T');
zlabel('\mu');
title('Plot of \mu(C, T)');