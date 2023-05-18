clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Simulation parameters
dh = 1e-3;
lambda_arr = [0.1, 0.2, 0.4, 0.6, 0.8];
plot_m = 3;
plot_n = 2;
% lambda_arr = [0.4];
% plot_m = 1;
% plot_n = 2;
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

p_gt_fun = test_case_data.p_gt_fun;
x_axis = linspace(left_first, right_first, N);
p_first = p_gt_fun(x_axis, 0)';

f = figure();
position = [100, 100, 500 * plot_n, 200 * plot_m];
set(f, 'Position', position);

subplot(plot_m, plot_n, 1);
hold on;
plot(x_axis, p_first);
xlim([pos_first-sigma, pos_first+sigma]);
ylim([0, 12]);
xlabel("x");
ylabel(sprintf("p(x,t=%.1f)", 0));
title(sprintf("Wave packet at t = %.1f", 0));


fft_size = 1024;
f_max = 1/dh;
f_axis = linspace(-f_max/2,f_max/2-dh,fft_size);

fft_first = fftshift(fft(p_first, fft_size));

f2 = figure();
set(f2, 'Position', position);
subplot(plot_m, plot_n, 1);
plot(f_axis, unwrap(angle(fft_first)));
xlim([-f_max/2,f_max/2-dh]);
ylim([-360,60]);
xlabel("f");
ylabel(sprintf("FFT{p(x,t=%.1f)}", 0));
title(sprintf("Wave packet at t = %.1f", 0));


pos_last = len_x/2 + c*1;
left_last = pos_last - sigma;
right_last = pos_last + sigma;

for n = 1:length(lambda_arr)
    lambda = lambda_arr(n);
    dt = dh * lambda / c;
    % Run simulation
    [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu, nu_fourier);

    p_last = p(round(left_last/dh):round(right_last/dh),end);

    str = sprintf("\\lambda = %.2f", lambda);
    
    figure(f);
    subplot(plot_m, plot_n, 1+n);
    plot(x_axis(round(left_last/dh):round(right_last/dh)), p_last, DisplayName=str);
    title(sprintf("Wave packet at t = %.1f - %s", len_t, str));
    xlim([pos_last-sigma, pos_last+sigma]);
    ylim([0, 12]);
    xlabel("x");
    ylabel(sprintf("p(x,t=%.1f)", len_t));

    fft_last = fftshift(fft(p_last, fft_size));

    figure(f2);
    subplot(plot_m, plot_n, 1+n);
    plot(f_axis, unwrap(angle(fft_last./fft_first)));
    title(sprintf("Wave packet at t = %.1f - %s", len_t, str));
    xlim([-f_max/2,f_max/2-dh]);
    ylim([-100,100]);
    xlabel("f");
    ylabel(sprintf("FFT{p(x,t=%.1f)}", len_t));
end








