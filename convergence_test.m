clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();

% Extracting test case data
p_gt_fun = test_case_data.p_gt_fun;
v_gt_fun = test_case_data.v_gt_fun;

% Simulation parameters
dt_values = [0.02 0.01 0.005 0.0025 0.001 0.0005];
dh_values = [0.2 0.1 0.05 0.025 0.01 0.005];

L1Err = zeros(length(dt_values), length(dh_values));
L2Err = zeros(length(dt_values), length(dh_values));
LinfErr = zeros(length(dt_values), length(dh_values));

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
stem(dt_values, L1Err);
set(gca, 'XScale', 'log');
title('L1 Error');
xlabel('dt');
ylabel('Error');

subplot(3,1,2);
stem(dt_values, L2Err);
set(gca, 'XScale', 'log');
title('L2 Error');
xlabel('dt');
ylabel('Error');

subplot(3,1,3);
stem(dt_values, LinfErr);
set(gca, 'XScale', 'log');
title('Linf Error');
xlabel('dt');
ylabel('Error');







