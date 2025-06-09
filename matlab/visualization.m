%----------------------------------------------------------------%
% __________        .__                __________                %
% \______   \_______|__|____    ____   \______   \ ____   ____   %
%  |    |  _/\_  __ \  \__  \  /    \   |       _// __ \ /    \  %
%  |    |   \ |  | \/  |/ __ \|   |  \  |    |   \  ___/|   |  \ %
%  |______  / |__|  |__(____  /___|  /  |____|_  /\___  >___|  / %
%         \/                \/     \/          \/     \/     \/  %
%                                                                %
%           Quasi-one-dimensional DeLaval Nozzle Flow Solver     %
%           Author : Zhixin Ren                                  %
%           Date   : 2025/5/24                                   %
%           License: MIT                                         %
%----------------------------------------------------------------%

close all;

fig_width = 30;
% Plot residuals
figure;
set(gcf,'Units', 'centimeters', 'OuterPosition',[1 1 fig_width 20]);
semilogy(1:max_iter, res(:, 1), 'LineWidth', 2);
hold on;
semilogy(1:max_iter, res(:, 2), 'LineWidth', 2);
hold on;
semilogy(1:max_iter, res(:, 3), 'LineWidth', 2);
xlabel('Iteration');
ylabel('Residual');
legend("Density", "Momentum", "Energy");
grid on;
set(gca,'OuterPosition',[.0 .0 1 1], 'FontSize',12, 'FontName', "NewComputerModern10");
saveas(gcf, 'residuals.svg');

% Plot density, temperature, and Mach number at the final time step
figure;
set(gcf,'Units', 'centimeters', 'OuterPosition',[fig_width+2, 22, fig_width 10]);
subplot(1, 3, 1);
plot(x, r(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Density");
set(gca,'OuterPosition',[.0 .0 .33 1], 'FontSize', 12, 'FontName', "NewComputerModern10");
grid on;
subplot(1, 3, 2);
plot(x, T(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Temperature");
set(gca,'OuterPosition',[.33 0 .33 1], 'FontSize', 12, 'FontName', "NewComputerModern10");
grid on;
subplot(1, 3, 3);
plot(x, M(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Mach number");
grid on;
set(gca,'OuterPosition',[.66 0 .33 1], 'FontSize', 12, 'FontName', "NewComputerModern10");
saveas(gcf, 'final_state.svg');

% Plot density, temperature, and Mach number at the final time step
figure;
set(gcf,'Units', 'centimeters', 'OuterPosition',[fig_width+2, 5, fig_width, fig_width]);
subplot(2, 2, 1);
plot(x, r(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Density");
set(gca,'OuterPosition',[.0 .5 .5 .5], 'FontSize', 12, 'FontName', "NewComputerModern10");
grid on;

subplot(2, 2, 2);
plot(x, T(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Temperature");
set(gca,'OuterPosition',[.5 .5 .5 .5], 'FontSize', 12, 'FontName', "NewComputerModern10");
grid on;

subplot(2, 2, 3);
plot(x, M(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Mach number");
grid on;

set(gca,'OuterPosition',[0 0 .5 .5], 'FontSize', 12, 'FontName', "NewComputerModern10");
subplot(2, 2, 4);
plot(x, p(end, :), 'LineWidth', 2);
xlabel('x');
ylabel("Pressure");
grid on;
set(gca,'OuterPosition',[.5 0 .5 .5], 'FontSize', 12, 'FontName', "NewComputerModern10");
saveas(gcf, 'final_state.svg');

% Plot mass flow rate through the nozzle
mass_flow = r .* V .* A;
figure;
set(gcf,'Units', 'centimeters', 'OuterPosition',[fig_width+2, 1, fig_width 20]);
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
legend("0dt", "50dt", "100dt", "150dt", "200dt", "700dt", "2000dt");
grid on;
set(gca,'OuterPosition',[.0 .0 1 1], 'FontSize', 12, 'FontName', "NewComputerModern10");
saveas(gcf, 'mass_flow_rate.svg');

% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');

