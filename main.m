clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();
c = test_case_data.c0;


% Plotting ground truth
plot_ground_truth(test_case_data);


% Simulation parameters
dh = 1e-1;
dt = dh * 0.05 / c;
% artificial dissipation factors (used for first order)
xi = 0.99;
nu = 0.9;

% Show debug plot?
debug = false;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, xi, nu);

% Plotting simulation
[fig_p, fig_v] = plot_spacetime(t_axis,x_axis,p,v,'Simulation');

% Save simulation as figures and animation
save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v);







