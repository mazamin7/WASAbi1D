clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Plotting ground truth
plot_ground_truth(test_case_data);


% Simulation parameters
dh = 1e-1;
dt = dh * 0.4 / c;

% artificial dissipation factors for first order
xi = 1 - 5e-2;
nu = 1; % 0.99; % in case of Fourier, it only affects DD

% Fourier artificial dissipation factor
nu_fourier = 1 - 1e-10; % - eps(1); % < 1

% Show debug plot?
debug = false;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu, nu_fourier);

% Plotting simulation
[fig_p, fig_v] = plot_spacetime(t_axis,x_axis,p,v,'Simulation');

% Save simulation as figures and animation
save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v);







