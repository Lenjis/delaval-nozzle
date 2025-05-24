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
plot(x, mass_flow(end, :), 'LineWidth', 2);
hold on;
xlabel('x');
ylabel('Mass Flow Rate');
title('Mass Flow Rate Along the Nozzle');
legend("0dt", "50dt", "100dt", "150dt", "200dt", "700dt", "final result");
grid on;
