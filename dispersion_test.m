clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Simulation parameters
dh = 1e-3;
% lambda_arr = [0.1, 0.2, 0.4, 0.6, 0.8];
% plot_m = 3;
% plot_n = 2;
lambda_arr = [1];
plot_m = 1;
plot_n = 2;
% lambda_arr = [0.6, 0.8];

% artificial dissipation factors for first order
xi = 1 - 1e-5;
nu = 1; % 0.99; % in case of Fourier, it only affects DD

% Fourier artificial dissipation factor
nu_fourier = 1 - 1e-10; % - eps(1); % < 1

% Show debug plot?
debug = false;

len_x = test_case_data.len_x;
len_t = test_case_data.len_t;
sigma = len_x/80;

pos_first = len_x/2;
left_first = pos_first - sigma;
right_first = pos_first + sigma;
N = round(right_first/dh) - round(left_first/dh) + 1;

pos_last = len_x/2 + c*1;
left_last = pos_last - sigma;
right_last = pos_last + sigma;

p_gt_fun = test_case_data.p_gt_fun;
N2 = 1e5;
x_axis2 = linspace(0, len_x, N2);
dh2 = len_x/N2;
p_gt = p_gt_fun(x_axis2, 0);
p_first = p_gt(round(left_first/dh2):round(right_first/dh2));

f = figure();
position = [100, 100, 500 * plot_n, 200 * plot_m];
set(f, 'Position', position);

subplot(plot_m, plot_n, 1);
hold on;
xlim([pos_first-sigma, pos_first+sigma]);
ylim([0, 12]);
xlabel("x");
ylabel(sprintf("p(x,t=%.1f)", 0));
title(sprintf("Wave packet at t = %.1f", 0));
plot(x_axis2(round(left_first/dh2):round(right_first/dh2)), p_first);

for n = 1:length(lambda_arr)
    lambda = lambda_arr(n);
    dt = dh * lambda / c;
    % Run simulation
    [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu, nu_fourier);

    p_last = p(round(left_last/dh):round(right_last/dh),end);

    str = sprintf("\\lambda = %.2f", lambda);
    
    subplot(plot_m, plot_n, 1+n);
    hold on;
    xlim([pos_last-sigma, pos_last+sigma]);
    ylim([0, 12]);
    xlabel("x");
    ylabel(sprintf("p(x,t=%.1f)", len_t));
    title(sprintf("Wave packet at t = %.1f - %s", len_t, str));
    plot(x_axis(round(left_last/dh):round(right_last/dh)), p_last, DisplayName=str);
end








