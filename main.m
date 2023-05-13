clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Plotting ground truth
plot_ground_truth(test_case_data);


% Simulation parameters
dh = 1e-1;
dt = dh * 0.3 / c;

% artificial dissipation factors for first order FDTD
xi = 0.999;
nu = 0.99;

% Fourier artificial dissipation factor
nu_fourier = 0.999999999; % < 1

% Show debug plot?
debug = false;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu, nu_fourier);

% Plotting simulation
[fig_p, fig_v] = plot_spacetime(t_axis,x_axis,p,v,'Simulation');

% Save simulation as figures and animation
save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v);







