function [test_case_data] = get_test_case()

% User menu
msg = "Choose the test case";
opts = ["1" "2" "3" "4" "5"];
test_case = menu(msg, opts);

switch test_case
    
    case 1 % gaussian shaped force
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

        % Defining force
        force_fun = @(x, t) force_time_fun(t) * force_envelope(x);
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0;
        g2_time_fun = @(t) 0;

        % Defining initial conditions
        p_gt_fun = @(x,t) 0 * x * t;
        v_gt_fun = @(x,t) 0 * x * t;
    
    case 2 % "string" driven on the left
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
        
        % Defining force
        force_fun = @(x, t) 0;
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0.1*cos(2*pi*1*t);
        g2_time_fun = @(t) 0;

        % Defining initial conditions
        p_gt_fun = @(x,t) 0 * x * t;
        v_gt_fun = @(x,t) 0 * x * t;

    case 3 % standing wave
        len_x = 1; % Domain length
        len_t = 2; % Simulation duration

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
        
        % Defining ground truth solution
        k_force = pi;
        omega_force = k_force * c0;
        p_gt_fun = @(x,t) 0.1*sin(omega_force*t) * cos(k_force * x);
        v_gt_fun = @(x,t) 0.1*omega_force*cos(omega_force*t) * cos(k_force * x);

        % Defining force
        force_fun = @(x, t) 0;
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0;
        g2_time_fun = @(t) 0;

    case 4 % propagating wave
        len_x = 10; % Domain length
        len_t = 4; % Simulation duration

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

        % Defining wave envelope
        source_sigma_ratio_x = 1/20;
        sigma = len_x * source_sigma_ratio_x;       % standard deviation of force spatial envelope (Gaussian)
        wave_envelope = @(x) 1/(sigma * sqrt(2 * pi)) * exp(-x.^2/(2*sigma^2)); % Gaussian function
        wave_envelope_d = @(x) -x/(sigma^2) .* wave_envelope(x);
        
        % Defining ground truth solution
        k_force = pi;
        omega_force = k_force * c0;
        p_gt_fun = @(x,t) 0.1 * ( wave_envelope(omega_force*t - k_force*(x-len_x/2)) + wave_envelope(omega_force*t + k_force*(x-len_x/2)) );
        v_gt_fun = @(x,t) 0.1 * omega_force * ( wave_envelope_d(omega_force*t - k_force*(x-len_x/2)) + wave_envelope_d(omega_force*t + k_force*(x-len_x/2)) );

        % Defining force
        force_fun = @(x, t) 0;
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0;
        g2_time_fun = @(t) 0;

    case 5 % propagating wave
        len_x = 10; % Domain length
        len_t = 4; % Simulation duration

        % Speed of propagation
        c0 = 1;

        % Absorption coefficients
        alpha_abs_left = 0.5;
        alpha_abs_right = 0.5;

        % Transmittance of the middle boundary
        transmittivity = 1;
        
        % Boundary conditions
        bc_left = "N";
        bc_right = "N";

        % Defining wave envelope
        source_sigma_ratio_x = 1/20;
        sigma = len_x * source_sigma_ratio_x;       % standard deviation of force spatial envelope (Gaussian)
        wave_envelope = @(x) 1/(sigma * sqrt(2 * pi)) * exp(-x.^2/(2*sigma^2)); % Gaussian function
        wave_envelope_d = @(x) -x/(sigma^2) .* wave_envelope(x);
        
        % Defining ground truth solution
        k_force = pi;
        omega_force = k_force * c0;
        p_gt_fun = @(x,t) 0.1 * exp(-alpha_abs_left*t) .* ... 
            ( wave_envelope(omega_force*t - k_force*(x-len_x/2)) ...
            + wave_envelope(omega_force*t + k_force*(x-len_x/2)) );
        v_gt_fun = @(x,t) 0.1 * exp(-alpha_abs_left*t) .* omega_force .* ...
            ( wave_envelope_d(omega_force*t - k_force*(x-len_x/2)) ...
            + wave_envelope_d(omega_force*t + k_force*(x-len_x/2)) );

        % Defining force
        force_fun = @(x, t) 0;
        
        % Defining boundary conditions time evolution
        g1_time_fun = @(t) 0;
        g2_time_fun = @(t) 0;

    otherwise
        error("Invalid test case number");
        
end

test_case_data.test_case = test_case;
test_case_data.len_x = len_x;
test_case_data.len_t = len_t;
test_case_data.c0 = c0;
test_case_data.alpha_abs_left = alpha_abs_left;
test_case_data.alpha_abs_right = alpha_abs_right;
test_case_data.transmittivity = transmittivity;
test_case_data.bc_left = bc_left;
test_case_data.bc_right = bc_right;
test_case_data.force_fun = force_fun;
test_case_data.g1_time_fun = g1_time_fun;
test_case_data.g2_time_fun = g2_time_fun;
test_case_data.p_gt_fun = p_gt_fun;
test_case_data.v_gt_fun = v_gt_fun;

end
