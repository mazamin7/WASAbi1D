clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();

% Extracting test case data
p_gt_fun = test_case_data.p_gt_fun;
v_gt_fun = test_case_data.v_gt_fun;

% Simulation parameters
dt_values = [0.02 0.01 0.005 0.0025 0.001 0.0005];
dh_values = [0.2 0.1 0.05 0.025 0.01 0.005];

L1Err = zeros(length(dt_values), 1);
L2Err = zeros(length(dt_values), 1);
LinfErr = zeros(length(dt_values), 1);

for i = 1:length(dt_values)
    dt = dt_values(i);
    dh = dh_values(i);

    % Run simulation
    [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh);

    % Evaluate the ground truth on the grid
    [X, T] = meshgrid(x_axis, t_axis);
    p_gt = p_gt_fun(X, T)';
    v_gt = v_gt_fun(X, T)';

    % Compute errors
    L1Err(i) = norm(p(:)-p_gt(:),1)/numel(p);
    L2Err(i) = sqrt(norm(p(:)-p_gt(:),2)^2/numel(p));
    LinfErr(i) = norm(p(:)-p_gt(:),inf);
end

% Plot errors
figure;
subplot(3,1,1);
stem(dt_values(end:-1:1), L1Err(end:-1:1));
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('L1 Error');
xlabel('dt');
ylabel('Error');
xticks(dt_values(end:-1:1));
xticklabels(arrayfun(@(x) sprintf('%.1e', x), dt_values(end:-1:1), 'UniformOutput', false));

subplot(3,1,2);
stem(dt_values(end:-1:1), L2Err(end:-1:1));
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('L2 Error');
xlabel('dt');
ylabel('Error');
xticks(dt_values(end:-1:1));
xticklabels(arrayfun(@(x) sprintf('%.1e', x), dt_values(end:-1:1), 'UniformOutput', false));

subplot(3,1,3);
stem(dt_values(end:-1:1), LinfErr(end:-1:1));
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('Linf Error');
xlabel('dt');
ylabel('Error');
xticks(dt_values(end:-1:1));
xticklabels(arrayfun(@(x) sprintf('%.1e', x), dt_values(end:-1:1), 'UniformOutput', false));








