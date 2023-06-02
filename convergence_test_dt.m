clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;
method = simulation_parameters.method_left;

% Extracting test case data
p_gt_fun = test_case_data.p_gt_fun;
v_gt_fun = test_case_data.v_gt_fun;

% Simulation parameters
dh = 1e-1;

if method == 1
    % FDTD 2ord
    lambda_arr = [0.001 0.002 0.004 0.008 0.01 0.02 0.04 0.08 0.1 0.2 0.4 0.8];
elseif method == 2
    % FDTD 1ord
    lambda_arr = [0.001 0.002 0.004 0.008 0.01 0.02 0.04 0.08 0.1 0.2 0.4 0.8]/2;
elseif method == 3
    % Fourier 2ord
    lambda_arr = [0.5 1 1.5 2];
elseif method == 4
    % Fourier 1ord
    lambda_arr = [0.5 1 1.5 2];
end

dt_arr = dh * lambda_arr / c;

% artificial dissipation factors for first order
xi = 1 - 1e-10;
nu = 1; % 0.99; % in case of Fourier, it only affects DD

% Fourier artificial dissipation factor
nu_fourier = 1 - 1e-10; % - eps(1); % < 1

L1Err = zeros(length(lambda_arr), 1);
L2Err = zeros(length(lambda_arr), 1);
LinfErr = zeros(length(lambda_arr), 1);

for i = 1:length(lambda_arr)
    dt = dt_arr(i);

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
plot(dt_arr, L1Err);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
hold on;

% Compute theoretical convergence curve
max_dt = max(dt_arr);
theory_curve = L1Err(end) * (dt_arr / max_dt);
plot(dt_arr, theory_curve, 'r-');
hold off;

title('L1 Error');
xlabel('dt');
ylabel('Error');
legend('Numerical', 'Theory');

subplot(3,1,2);
plot(dt_arr, L2Err);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
hold on;

% Compute theoretical convergence curve
theory_curve = L2Err(end) * (dt_arr / max_dt);
plot(dt_arr, theory_curve, 'r-');
hold off;

title('L2 Error');
xlabel('dt');
ylabel('Error');
legend('Numerical', 'Theory');

subplot(3,1,3);
plot(dt_arr, LinfErr);
set(gca, 'XScale', 'log', 'XDir', 'reverse');
hold on;

% Compute theoretical convergence curve
theory_curve = LinfErr(end) * (dt_arr / max_dt);
plot(dt_arr, theory_curve, 'r-');
hold off;

title('Linf Error');
xlabel('dt');
ylabel('Error');
legend('Numerical', 'Theory');









