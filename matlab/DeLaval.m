%----------------------------------------------------------------%
% __________        .__                __________                %
% \______   \_______|__|____    ____   \______   \ ____   ____   %
%  |    |  _/\_  __ \  \__  \  /    \   |       _// __ \ /    \  %
%  |    |   \ |  | \/  |/ __ \|   |  \  |    |   \  ___/|   |  \ %
%  |______  / |__|  |__(____  /___|  /  |____|_  /\___  >___|  / %
%         \/                \/     \/          \/     \/     \/  %
%                                                                %
%           Quasi-one-dimensional DeLaval Nozzle Flow Solver     %
%               Author : Zhixin Ren                              %
%               Date   : 2025/5/24                               %
%               License: MIT                                     %
%----------------------------------------------------------------%
clear;
close all;
C = 0.5; % Courant number
max_iter = 1000; % Maximum number of iterations
gamma = 1.4; % Specific heat ratio
% Geometry
dx = 0.1;
x = 0:dx:3;
A = 1 + 2.2 * (x - 1.5) .^ 2;

% results
r = zeros(max_iter + 1, length(x));
T = zeros(max_iter + 1, length(x));
V = zeros(max_iter + 1, length(x));
M = zeros(max_iter + 1, length(x));
% predicted values
r_p = zeros(max_iter, length(x));
T_p = zeros(max_iter, length(x));
V_p = zeros(max_iter, length(x));
% predicted time derivatives
r_t_p = zeros(max_iter, length(x));
T_t_p = zeros(max_iter, length(x));
V_t_p = zeros(max_iter, length(x));
% corrected time derivatives
r_t_c = zeros(max_iter, length(x));
T_t_c = zeros(max_iter, length(x));
V_t_c = zeros(max_iter, length(x));
% averaged time derivatives
r_t_a = zeros(max_iter, length(x));
T_t_a = zeros(max_iter, length(x));
V_t_a = zeros(max_iter, length(x));

res = zeros(max_iter, 3); % Residuals

% Initialization
r(1, :) = 1 - 0.3146 * x;
T(1, :) = 1 - 0.2314 * x;
V(1, :) = (0.1 + 1.09 * x) .* T(1, :) .^ 0.5;

% Boundary conditions
% Inlet: r and T are stipulated, V to float
% r(1, 1) = 1;
% T(1, 1) = 1;
% Outlet: All 3 variables are floated

% Solving using time-marching, MacCormack's explicit method
for t = 1:max_iter
    % Calculating dt
    dts = C * dx ./ (T(t, 2:end - 1) .^ 0.5 + V(t, 2:end - 1));
    dt = min(dts);

    r(t + 1, 1) = r(t, 1);
    T(t + 1, 1) = T(t, 1);

    % prediction step, forward difference in space
    for i = 1:(length(x) - 1)
        r_t_p(t, i) = -r(t, i) * (V(t, i + 1) - V(t, i)) / dx - r(t, i) * V(t, i) * (log(A(i + 1)) - log(A(i))) / dx - V(t, i) * (r(t, i + 1) - r(t, i)) / dx;
        V_t_p(t, i) = -V(t, i) * (V(t, i + 1) - V(t, i)) / dx - ((T(t, i + 1) - T(t, i)) / dx + T(t, i) / r(t, i) * (r(t, i + 1) - r(t, i)) / dx) / gamma;
        T_t_p(t, i) = -V(t, i) * (T(t, i + 1) - T(t, i)) / dx - (gamma - 1) * T(t, i) * ((V(t, i + 1) - V(t, i)) / dx + V(t, i) * (log(A(i + 1)) - log(A(i))) / dx);

        r_p(t, i) = r(t, i) + r_t_p(t, i) * dt;
        V_p(t, i) = V(t, i) + V_t_p(t, i) * dt;
        T_p(t, i) = T(t, i) + T_t_p(t, i) * dt;
    end

    r_p(t, 1) = r(t, 1);
    T_p(t, 1) = T_p(t, 1);

    % correction step, rearward difference in space
    for i = 2:(length(x) - 1)
        r_t_c(t, i) = -r_p(t, i) * (V_p(t, i) - V_p(t, i - 1)) / dx - r_p(t, i) * V_p(t, i) * (log(A(i)) - log(A(i - 1))) / dx - V_p(t, i) * (r_p(t, i) - r_p(t, i - 1)) / dx;
        V_t_c(t, i) = -V_p(t, i) * (V_p(t, i) - V_p(t, i - 1)) / dx - ((T_p(t, i) - T_p(t, i - 1)) / dx + T_p(t, i) / r_p(t, i) * (r_p(t, i) - r_p(t, i - 1)) / dx) / gamma;
        T_t_c(t, i) = -V_p(t, i) * (T_p(t, i) - T_p(t, i - 1)) / dx - (gamma - 1) * T_p(t, i) * ((V_p(t, i) - V_p(t, i - 1)) / dx + V_p(t, i) * (log(A(i)) - log(A(i - 1))) / dx);

        r_t_a(t, i) = (r_t_p(t, i) + r_t_c(t, i)) / 2;
        V_t_a(t, i) = (V_t_p(t, i) + V_t_c(t, i)) / 2;
        T_t_a(t, i) = (T_t_p(t, i) + T_t_c(t, i)) / 2;

        r(t + 1, i) = r(t, i) + r_t_a(t, i) * dt;
        V(t + 1, i) = V(t, i) + V_t_a(t, i) * dt;
        T(t + 1, i) = T(t, i) + T_t_a(t, i) * dt;
    end

    M(t + 1, :) = V(t + 1, :) ./ T(t + 1, :) .^ 0.5; % Mach number

    % Linear extrapolation for boundary nodes
    V(t + 1, 1) = 2 * V(t + 1, 2) - V(t + 1, 3);
    r(t + 1, end) = 2 * r(t + 1, end - 1) - r(t + 1, end - 2);
    T(t + 1, end) = 2 * T(t + 1, end - 1) - T(t + 1, end - 2);
    V(t + 1, end) = 2 * V(t + 1, end - 1) - V(t + 1, end - 2);

    % Compute and plot residuals (L2 norm of change in variables)
    res(t, :) = [norm(r(t + 1, :) - r(t, :)), norm(T(t + 1, :) - T(t, :)), norm(V(t + 1, :) - V(t, :))];

end

% Visualization
figure;
semilogy(1:max_iter, res(:, 1), 'LineWidth', 2);
hold on;
semilogy(1:max_iter, res(:, 2), 'LineWidth', 2);
hold on;
semilogy(1:max_iter, res(:, 3), 'LineWidth', 2);
xlabel('Iteration');
title('Residual');

figure;
subplot(1, 3, 1);
plot(x, r(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Density");
subplot(1, 3, 2);
plot(x, T(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Temperature");
subplot(1, 3, 3);
plot(x, M(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Mach number");

figure;
subplot(1, 3, 1);
plot(1:max_iter, r(2:end, 16), 'LineWidth', 2);
xlabel('Iteration');
ylabel("Density");
subplot(1, 3, 2);
plot(1:max_iter, T(2:end, 16), 'LineWidth', 2);
xlabel('Iteration');
ylabel("Temperature");
subplot(1, 3, 3);
plot(1:max_iter, M(2:end, 16), 'LineWidth', 2);
xlabel('Iteration');
ylabel("Mach number");

% Plot mass flow rate through the nozzle
mass_flow = r .* V .* A;
figure;
plot(x, mass_flow(1, :), 'LineWidth', 2);
hold on;
plot(x, mass_flow(51, :), 'LineWidth', 2);
hold on;
plot(x, mass_flow(101, :), 'LineWidth', 2);
hold on;
plot(x, mass_flow(151, :), 'LineWidth', 2);
hold on;
plot(x, mass_flow(201, :), 'LineWidth', 2);
hold on;
plot(x, mass_flow(701, :), 'LineWidth', 2);
hold on;
xlabel('x');
ylabel('Mass Flow Rate');
title('Mass Flow Rate Along the Nozzle');
legend("0dt", "50dt", "100dt", "150dt", "200dt", "700dt");
grid on;
