clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Simulation parameters
dh = 1e-3;
dt = dh * 0.4 / c;

% artificial dissipation factors for first order
xi = 1 - 1e-5;
nu = 1; % 0.99; % in case of Fourier, it only affects DD

% Fourier artificial dissipation factor
nu_fourier = 1 - 1e-10; % - eps(1); % < 1

% Show debug plot?
debug = false;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu, nu_fourier);

len_x = test_case_data.len_x;
len_t = test_case_data.len_t;
sigma = 3 * 1/20 * len_x;

figure();
sgtitle("Wave packet dispersion test");

subplot(211);
pos_first = len_x/2;
left_first = pos_first - sigma;
right_first = pos_first + sigma;
p_first = p(round(left_first/dh):round(right_first/dh),1);
N = round(right_first/dh) - round(left_first/dh) + 1;
plot(x_axis(round(left_first/dh):round(right_first/dh)), p_first);
xlim([pos_first-len_x/40, pos_first+len_x/40]);
ylim([0, 12]);
xlabel("x");
ylabel(sprintf("p(x,t=%.1f)", 0));
title(sprintf("Wave packet at t = %.1f", 0));

subplot(212);
pos_last = len_x/2 + c*1;
left_last = pos_last - sigma;
right_last = pos_last + sigma;
p_last = p(round(left_last/dh):round(right_last/dh),end);
plot(x_axis(round(left_last/dh):round(right_last/dh)), p_last);
xlim([pos_last-len_x/40, pos_last+len_x/40]);
ylim([0, 12]);
xlabel("x");
ylabel(sprintf("p(x,t=%.1f)", len_t));
title(sprintf("Wave packet at t = %.1f", len_t));

% fft_len = 8192;
% % Compute the Fourier transform of the first and last time instant
% fft_first = fft(p_first,fft_len);
% fft_last = fft(p_last,fft_len);
% 
% % Define the frequency axis
% f_axis = -fft_len/2:fft_len/2-1;  % Frequency axis
% 
% figure;
% subplot(2, 1, 1);
% plot(f_axis, abs(fft_first));
% subplot(2, 1, 2);
% plot(f_axis, abs(fft_last));

% figure;
% subplot(2, 1, 1);
% plot(f_axis, ifft(fft_first));
% subplot(2, 1, 2);
% plot(f_axis, ifft(fft_last));

% % Compute the impulse response using the inverse Fourier transform
% freq_resp = fft_last ./ fft_first;
% 
% % Plotting the frequency response
% figure;
% plot(f_axis, abs(freq_resp));
% title('Frequency Response');
% xlabel('Frequency');
% ylabel('Magnitude');








