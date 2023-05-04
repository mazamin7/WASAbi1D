clear all, close all, clc;

simulation_parameters = get_simulation_parameters();
test_case_data = get_test_case();


% Plotting ground truth
plot_ground_truth(test_case_data);


% Simulation parameters
dt = 0.005;
dh = 0.05;

% Run simulation
[t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, true);

% Plotting simulation
[fig_p, fig_v] = plot_spacetime(t_axis,x_axis,p,v,'Simulation');

% Save simulation as figures and animation
save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v);






