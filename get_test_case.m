function [len_x, len_t, c0, alpha_abs_left, alpha_abs_right, transmittivity, bc_left, bc_right, force_time_fun, force_envelope, g1_time_fun, g2_time_fun] = get_test_case(test_case)

switch test_case
    
    case 1
        len_x = 10; % Domain length
        len_t = 10; % Simulation duration

        % Speed of propagation
        c0 = 1;

        % Absorption coefficients
        alpha_abs_left = 0;
        alpha_abs_right = 0;

        % Transmittance of the middle boundary
        transmittivity = 1;
        
        % Boundary conditions
        bc_left = "N";
        bc_right = "N";
        
        % Defining force time evolution
        freq_force = 1;
        force_time_fun = @(t) sin(2*pi*freq_force*t) * (t <= 1/freq_force);
        
        % Defining force spatial envelope
        source_mu_ratio_x = 1/4;
        source_sigma_ratio_x = 1/40;
        mu = len_x * source_mu_ratio_x;
        sigma = len_x * source_sigma_ratio_x;       % standard deviation of force spatial envelope (Gaussian)
        force_envelope = @(x) 1/(sigma * sqrt(2 * pi)) * exp(-(x-mu).^2/(2*sigma^2)); % Gaussian function
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0;
        g2_time_fun = @(t) 0;
    
    case 2
        len_x = 10; % Domain length
        len_t = 10; % Simulation duration

        % Speed of propagation
        c0 = 1;

        % Absorption coefficients
        alpha_abs_left = 0.3;
        alpha_abs_right = 1;

        % Transmittance of the middle boundary
        transmittivity = 1;
        
        % Boundary conditions
        bc_left = "D";
        bc_right = "D";
        
        % Defining force time evolution
        force_time_fun = @(t) 0;
        
        % Defining force spatial envelope
        force_envelope = @(x) 0; % Gaussian function
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0.1*cos(2*pi*1*t);
        g2_time_fun = @(t) 0;

    otherwise
        error("Invalid test case number.");
        
end
