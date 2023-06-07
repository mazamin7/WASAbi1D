function [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, diss, db_plot)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here

    % Extracting simulation parameters
    merge = simulation_parameters.merge;
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    order_left = simulation_parameters.order_left;
    order_right = simulation_parameters.order_right;
    DD = simulation_parameters.DD;

    % Extracting test case data
    len_x = test_case_data.len_x;
    len_t = test_case_data.len_t;
    c0 = test_case_data.c0;
    alpha_abs = test_case_data.alpha_abs;
    transmittivity = test_case_data.transmittivity;
    bc_left = test_case_data.bc_left;
    bc_right = test_case_data.bc_right;
    force_fun = test_case_data.force_fun;
    g1_time_fun = test_case_data.g1_time_fun;
    g2_time_fun = test_case_data.g2_time_fun;
    p_gt_fun = test_case_data.p_gt_fun;
    v_gt_fun = test_case_data.v_gt_fun;


    fourier_left = method_left == 3 || method_left == 4;
    fourier_right = method_right == 3 || method_right == 4;

    
    % Checking compatibility between simulation parameters and test case
    assert(~((method_left == 3 || method_left == 4) && (bc_left == "D" || bc_right == "D")), ...
        'Boundary conditions must be NEUMANN homogeneous with Fourier method')
    
    damped = alpha_abs ~= 0;
    
    assert(~((method_left == 3 || method_right == 3) && damped), 'Fourier 2ord does not support damping');
    
    if DD == true
        stable = check_stability(c0, dt, dh, order_left, fourier_left, diss, DD);
        assert(stable, 'Stability condition not satisfied on the left');
    
        stable = check_stability(c0, dt, dh, order_right, fourier_right, diss, DD);
        assert(stable, 'Stability condition not satisfied on the right');
    elseif DD == false
        stable = check_stability(c0, dt, dh, order_left, fourier_left, diss, DD);
        assert(stable, 'Stability condition not satisfied');
    end
    
    % Knowing simulation pars and test case, initialize simulation variables
    
    % Defining time and space axes
    N_t = floor(len_t / dt);
    N_x = floor(len_x / dh);
    
    x_axis = linspace(0,len_x,N_x);
    t_axis = linspace(0,len_t,N_t);
    
    % Initialize solution, force, and boundary conditions data
    p = zeros(N_x,N_t);
    v = zeros(N_x,N_t);
    force = zeros(N_x,N_t);
    g1 = zeros(N_t,1);
    g2 = zeros(N_t,1);
    p_gt = zeros(N_x,N_t);
    v_gt = zeros(N_x,N_t);
    
    for n = 1:N_t
        force(:,n) = force_fun(x_axis, t_axis(n));
        g1(n) = g1_time_fun(t_axis(n));
        g2(n) = g2_time_fun(t_axis(n));
        p_gt(:,n) = p_gt_fun(x_axis, t_axis(n));
        v_gt(:,n) = v_gt_fun(x_axis, t_axis(n));
    end
    
    % Imposing initial conditions
    p(:,1:2) = p_gt(:,1:2);
    v(:,1:2) = v_gt(:,1:2);
    
    if (method_left == 3 || method_left == 4 || method_right == 3 || method_right == 4)
        assert(all(g1 == 0), 'Boundary conditions must be Neumann HOMOGENEOUS with Fourier method');
        assert(all(g2 == 0), 'Boundary conditions must be Neumann HOMOGENEOUS with Fourier method');
    end

    if DD
        % Building residue matrix
        C = get_residue_matrix(N_x, 6);
        
        % Initializing update methods
        if fourier_left == false
            data_left = init_FDTD(len_x/2, c0, dt, dh, alpha_abs, bc_left, "N", method_left == 5, order_left);
        else
            data_left = init_Fourier(len_x/2, c0, dt, dh, order_left, alpha_abs);
        end
        
        if fourier_right == false
            data_right = init_FDTD(len_x/2, c0, dt, dh, alpha_abs, "N", bc_right, method_right == 5, order_right);
        else
            data_right = init_Fourier(len_x/2, c0, dt, dh, order_right, alpha_abs);
        end
    else
        % Initializing update methods
        if fourier_left == false
            data_left = init_FDTD(len_x, c0, dt, dh, alpha_abs, bc_left, bc_right, method_left == 5, order_left);
        else
            data_left = init_Fourier(len_x, c0, dt, dh, order_left, alpha_abs);
        end
    end

    if debug
        % Init figure
        f = figure();
        f.Position = [100, 100, 1200, 700];
    end

    force_n_offset = order_left == 1;% && fourier_left == false;
    
    % Simulation loop
    for n = 2:N_t-1
        % Step 0 (prediction correction, only if DD and pre-merge and first order)
%         if DD && order_left == 1 && merge == 1
%             % first guess
%             residual = (c0 * dt / dh)^2 * C * v(:,n) + (c0 / dh)^2 * C * p(:,n);
%             force_now = force(:,n+force_n_offset) + transmittivity^2 * residual;
% 
%             % iterative refinement
%             for it = 1:10 % min 1
%                 % Update left
%                 if method_left <= 2 || method_left == 5
%                     [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_FDTD(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n), g1(n), 0);
%                 elseif method_left >= 3
%                     [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_Fourier(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n));
%                 end
%     
%                 % Update right
%                 if method_right <= 2 || method_right == 5
%                     [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_FDTD(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n), 0, g2(n));
%                 elseif method_right >= 3
%                     [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_Fourier(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n));
%                 end
% 
%                 residual = (c0 / dh)^2 * C * p(:,n+1);
%                 force_now = force(:,n+force_n_offset) + transmittivity^2 * residual;
%             end
%         end

        if DD
            % Compute residual
            if order_left == 2
                residual = (c0 / dh)^2 * C * p(:,n); % exact
            elseif order_left == 1 && fourier_left == true
                % make a guess of p(n+1) using FDTD
                p(:,n+1) = dt * v(:,n) + p(:,n);
                residual = (c0 / dh)^2 * C * p(:,n+1);
            else % order_left == 1 && fourier_left == false
                % make a guess of p(n+1) using p(n)
                residual = (c0 / dh)^2 * C * p(:,n);
            end

            if order_left == 2 % Second order

                % Pre-merge
                if merge == 1
                    force_now = force(:,n+force_n_offset) + transmittivity^2 * residual;
                else
                    force_now = force(:,n+force_n_offset);
                end
    
                % Update left
                if method_left <= 2 || method_left == 5
                    [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_FDTD(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n), g1(n), 0);
                elseif method_left >= 3
                    [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_Fourier(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n));
                end
        
                % Update right
                if method_right <= 2 || method_right == 5
                    [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_FDTD(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n), 0, g2(n));
                elseif method_right >= 3
                    [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_Fourier(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n));
                end
                
                % Post-merge
                if merge == 2
                    p(:,n+1) = p(:,n+1) + transmittivity^2 * dt*dt * residual / (1 + dt*alpha_abs);
                end

            elseif order_left == 1 && merge == 2 % First order post merge

                force_now = force(:,n+force_n_offset);

                % Update left
                if method_left <= 2 || method_left == 5
                    [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_FDTD(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n), g1(n), 0);
                elseif method_left >= 3
                    [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_Fourier(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n));
                end
        
                % Update right
                if method_right <= 2 || method_right == 5
                    [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_FDTD(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n), 0, g2(n));
                elseif method_right >= 3
                    [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_Fourier(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n));
                end
                
                % Post-merge
                residual = (c0 / dh)^2 * C * p(:,n+1);
                v(:,n+1) = v(:,n+1) + transmittivity^2 * dt * residual / (1 + 2*dt*alpha_abs);

            elseif order_left == 1 && merge == 1 % First order pre merge

                % Prediction correction scheme

                % Decide how many iterations
                if fourier_left
                    num_it = 10;
                else
                    num_it = 2;
                end

                for it = 1:num_it
                    % Pre-merge
                    if merge == 1
                        force_now = force(:,n+force_n_offset) + transmittivity^2 * residual;
                    else
                        force_now = force(:,n+force_n_offset);
                    end
    
                    % Update left
                    if method_left <= 2 || method_left == 5
                        [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_FDTD(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n), g1(n), 0);
                    elseif method_left >= 3
                        [p(1:N_x/2,n+1),v(1:N_x/2,n+1)] = update_Fourier(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_now(1:N_x/2), v(1:N_x/2,n));
                    end
        
                    % Update right
                    if method_right <= 2 || method_right == 5
                        [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_FDTD(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n), 0, g2(n));
                    elseif method_right >= 3
                        [p(N_x/2+1:N_x,n+1),v(N_x/2+1:N_x,n+1)] = update_Fourier(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n));
                    end

                    % Recompute residual
                    residual = (c0 / dh)^2 * C * p(:,n+1);
                end

            end
        else % no DD
            force_now = force(:,n+force_n_offset);

            % Update
            if method_left <= 2 || method_left == 5
                [p(:,n+1),v(:,n+1)] = update_FDTD(data_left, p(:,n), p(:,n-1), force_now(:), v(:,n), g1(n), 0);
            elseif method_left >= 3
                [p(:,n+1),v(:,n+1)] = update_Fourier(data_left, p(:,n), p(:,n-1), force_now(:), v(:,n));
            end
        end

        % Artificial dissipation for stability
        p(:,n+1) = diss * p(:,n+1);
    
        % Computing velocity if second order
        if order_left == 2
            v(1:N_x/2,n+1) = (p(1:N_x/2,n+1) - p(1:N_x/2,n-1))/2/dt;
        end

        if order_right == 2
            v(N_x/2+1:N_x,n+1) = (p(N_x/2+1:N_x,n+1) - p(N_x/2+1:N_x,n-1))/2/dt;
        end
        
        if debug == true
            % Plot
            figure(f);
            sgtitle(['instant [s]: ' num2str((n+1)*dt, '%4.3f') ' / ' ...
                num2str(len_t, '%4.3f') ' ( ' num2str((n+1)/N_t*100, '%4.1f') '% )']);
        
            if mod(n,10) == 1
	            if db_plot == false
		            % Plot p
		            subplot(2,1,1);
		            plot(x_axis, p(:,n+1));
		            title('Pressure');
		            xlim([0,len_x]);
		            ylim([-1,1]*2e-1);
            
		            % Plot v
		            subplot(2,1,2);
		            plot(x_axis, v(:,n+1));
		            title('Velocity');
		            xlim([0,len_x]);
		            ylim([-c0,c0]*5e-1);
	            else
                    % Plot p
                    subplot(211);
                    plot(x_axis, db(p(:,n+1)));
                    title('Pressure');
                    grid on;
                    xlim([0,len_x]);
                    ylim([-150 0]);
                    yticks(-150:10:0);
            
                    % Plot v
                    subplot(212);
                    plot(x_axis, db(v(:,n+1)));
                    title('Velocity');
                    grid on;
                    xlim([0,len_x]);
                    ylim([-150 0]);
                    yticks(-150:10:0);
                end
            end
        else
            clc;
            disp(['Simulation: ' num2str((n+1)/N_t*100) '%']);
        end
    
    end
    
end

