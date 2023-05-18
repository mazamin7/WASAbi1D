clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;
method = simulation_parameters.method_left;

% Extracting test case data
p_gt_fun = test_case_data.p_gt_fun;
v_gt_fun = test_case_data.v_gt_fun;

% Simulation parameters
dh = 1e-2;

if method == 1
    % FDTD 2ord
    lambda_arr = [0.1 0.2 0.4 0.6 0.8];
elseif method == 2
    % FDTD 1ord
    lambda_arr = [0.05 0.1 0.15 0.2 0.25 0.3];
elseif method == 3
    % Fourier 2ord
    lambda_arr = [1];
elseif method == 4
    % Fourier 1ord
    lambda_arr = [0.5];
end

% artificial dissipation factors for first order
xi = 1 - 1e-10;
nu = 1; % 0.99; % in case of Fourier, it only affects DD

% Fourier artificial dissipation factor
nu_fourier = 1 - 1e-10; % - eps(1); % < 1

L1Err = zeros(length(lambda_arr), 1);
L2Err = zeros(length(lambda_arr), 1);
LinfErr = zeros(length(lambda_arr), 1);

for i = 1:length(lambda_arr)
    lambda = lambda_arr(i);
    dt = dh * lambda / c;

    % Run simulation
    [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, false, xi, nu, nu_fourier);

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
stem(lambda_arr, L1Err);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('L1 Error');
xlabel('dt');
ylabel('Error');
xticks(lambda_arr);
xticklabels(arrayfun(@(x) sprintf('%.1e', x), lambda_arr, 'UniformOutput', false));

subplot(3,1,2);
stem(lambda_arr, L2Err);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('L2 Error');
xlabel('dt');
ylabel('Error');
xticks(lambda_arr);
xticklabels(arrayfun(@(x) sprintf('%.1e', x), lambda_arr, 'UniformOutput', false));

subplot(3,1,3);
stem(lambda_arr, LinfErr);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('Linf Error');
xlabel('dt');
ylabel('Error');
xticks(lambda_arr);
xticklabels(arrayfun(@(x) sprintf('%.1e', x), lambda_arr, 'UniformOutput', false));








