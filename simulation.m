function [t_axis, x_axis, p, v] = simulation(test_case_data, simulation_parameters, dt, dh, debug, diss, db_plot)
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here

    % Extracting simulation parameters
    merge_left = simulation_parameters.merge_left;
    merge_right = simulation_parameters.merge_right;
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    space_order = simulation_parameters.space_order;
    temp_order_left = simulation_parameters.temp_order_left;
    temp_order_right = simulation_parameters.temp_order_right;
    DD = simulation_parameters.DD;

    % Extracting test case data
    len_x = test_case_data.len_x;
    len_t = test_case_data.len_t;
    c0 = test_case_data.c0;
    alpha_abs = test_case_data.alpha_abs;
    transmittivity = test_case_data.transmittivity;
    force_fun = test_case_data.force_fun;
    p_gt_fun = test_case_data.p_gt_fun;
    v_gt_fun = test_case_data.v_gt_fun;


    fourier_left = method_left == 3 || method_left == 4;
    fourier_right = method_right == 3 || method_right == 4;
    
    damped = alpha_abs ~= 0;
    
    assert(~((method_left == 3 || method_right == 3) && damped), 'Fourier 2ord does not support damping');
    
    if DD == true
        stable = check_stability(c0, dt, dh, temp_order_left, fourier_left, diss, DD);
        assert(stable, 'Stability condition not satisfied on the left');
    
        stable = check_stability(c0, dt, dh, temp_order_right, fourier_right, diss, DD);
        assert(stable, 'Stability condition not satisfied on the right');
    elseif DD == false
        stable = check_stability(c0, dt, dh, temp_order_left, fourier_left, diss, DD);
        assert(stable, 'Stability condition not satisfied');
    end
    
    % Knowing simulation pars and test case, initialize simulation variables
    
    % Defining time and space axes
    N_t = floor(len_t / dt);
    N_x = floor(len_x / dh);
    
    x_axis = linspace(0,len_x,N_x);
    t_axis = linspace(0,len_t,N_t);
    
    % Initialize solution, force, and boundary conditions data
    p = zeros(N_x,N_t+1);
    v = zeros(N_x,N_t+1);
    force = zeros(N_x,N_t+1);
    p_gt = zeros(N_x,N_t+1);
    v_gt = zeros(N_x,N_t+1);
    
    for n = 2:N_t+1
        force(:,n) = force_fun(x_axis, t_axis(n-1));
        p_gt(:,n) = p_gt_fun(x_axis, t_axis(n-1));
        v_gt(:,n) = v_gt_fun(x_axis, t_axis(n-1));
    end
    
    % Step 0
    % Imposing initial conditions (in position 1 we have stub)
    p(:,2) = p_gt(:,2);
    v(:,2) = v_gt(:,2);

    if DD
        % Building residue matrix
        C = get_residue_matrix(N_x, space_order);
        
        % Initializing update methods
        if fourier_left == false
            data_left = init_FDTD(len_x/2, c0, dt, dh, space_order, temp_order_left, alpha_abs);
        else
            data_left = init_Fourier(len_x/2, c0, dt, dh, temp_order_left, alpha_abs);
        end
        
        if fourier_right == false
            data_right = init_FDTD(len_x/2, c0, dt, dh, space_order, temp_order_right, alpha_abs);
        else
            data_right = init_Fourier(len_x/2, c0, dt, dh, temp_order_right, alpha_abs);
        end
    else
        % Initializing update methods
        if fourier_left == false
            data_left = init_FDTD(len_x, c0, dt, dh, space_order, temp_order_left, alpha_abs);
        else
            data_left = init_Fourier(len_x, c0, dt, dh, temp_order_left, alpha_abs);
        end
    end

    if debug
        % Init figure
        f = figure();
        f.Position = [100, 100, 1200, 700];
    end


    % Impose force
    force_next = force(:,2);

    if DD
        % Compute residual
        residual = (c0 / dh)^2 * C * p(:,2);
        
        % Pre-merge
        % Will be used in the next step for second order:
        %                    using r^n, on f^n
        %              in the current step for first order:
        %                    using r^{n+1}, on f^{n+1}
        force_next_corr = force_next * 0;

        if merge_left == 1
            force_next_corr(1:N_x/2) = force_next(1:N_x/2) + transmittivity^2 * residual(1:N_x/2);
        end

        if merge_right == 1
            force_next_corr(N_x/2+1:N_x) = force_next(N_x/2+1:N_x) + transmittivity^2 * residual(N_x/2+1:N_x);
        end
    end

    % Steps 1:N_t-1 (shifted by 1 to simplify code)
    
    % Simulation loop
    for n = 2:N_t
        override_order = n == 2;

        if DD

            % Update force
            force_now = force_next;
            force_now_corr = force_next_corr;


            % If order == 2, use pre-merge force in pressure
            force_use = force_now * 0;

            if temp_order_left == 1 || override_order
                force_use(1:N_x/2) = force_now(1:N_x/2);
            else
                force_use(1:N_x/2) = force_now_corr(1:N_x/2);
            end

            if temp_order_right == 1 || override_order
                force_use(N_x/2+1:N_x) = force_now(N_x/2+1:N_x);
            else
                force_use(N_x/2+1:N_x) = force_now_corr(N_x/2+1:N_x);
            end

            % Update pressure left
            if method_left <= 2 || method_left == 5
                % current force corrected in second order
                % (unused) current force not corrected in first order
                p(1:N_x/2,n+1) = update_pressure_FDTD(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_use(1:N_x/2), v(1:N_x/2,n), override_order);
            elseif method_left >= 3
                % current force corrected in second order
                % current force not corrected in first order
                p(1:N_x/2,n+1) = update_pressure_Fourier(data_left, p(1:N_x/2,n), p(1:N_x/2,n-1), force_use(1:N_x/2), v(1:N_x/2,n), override_order);
            end
            
            % Update pressure right
            if method_right <= 2 || method_right == 5
                p(N_x/2+1:N_x,n+1) = update_pressure_FDTD(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_use(N_x/2+1:N_x), v(N_x/2+1:N_x,n), override_order);
            elseif method_right >= 3
                p(N_x/2+1:N_x,n+1) = update_pressure_Fourier(data_right, p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_use(N_x/2+1:N_x), v(N_x/2+1:N_x,n), override_order);
            end
            
            % Post-merge second order
            % using r^n, on p^{n+1}
            if merge_left == 2
                if temp_order_left == 2 && override_order == false
                    p(1:N_x/2,n+1) = p(1:N_x/2,n+1) + transmittivity^2 * dt*dt * residual(1:N_x/2) / (1 + dt*alpha_abs);
                end
            end

            if merge_right == 2
                if temp_order_right == 2 && override_order == false
                    p(N_x/2+1:N_x,n+1) = p(N_x/2+1:N_x,n+1) + transmittivity^2 * dt*dt * residual(N_x/2+1:N_x) / (1 + dt*alpha_abs);
                end
            end
            
            % Compute new residual
            residual = (c0 / dh)^2 * C * p(:,n+1);
            
            % Impose new force
            force_next = force(:,n+1);

            % Artificial dissipation for stability
            p(:,n+1) = diss * p(:,n+1);
            

            % Pre-merge
            % Will be used in the next step for second order:
            %                    using r^n, on f^n
            %              in the current step for first order:
            %                    using r^{n+1}, on f^{n+1}
            if merge_left == 1
                force_next_corr(1:N_x/2) = force_next(1:N_x/2) + transmittivity^2 * residual(1:N_x/2);
            end

            if merge_right == 1
                force_next_corr(N_x/2+1:N_x) = force_next(N_x/2+1:N_x) + transmittivity^2 * residual(N_x/2+1:N_x);
            end


            % If order == 1, use pre-merge force in velocity
            force_use = force_now * 0;

            if temp_order_left == 2
                force_use(1:N_x/2) = force_now(1:N_x/2);
            else
                force_use(1:N_x/2) = force_next_corr(1:N_x/2);
            end

            if temp_order_right == 2
                force_use(N_x/2+1:N_x) = force_now(N_x/2+1:N_x);
            else
                force_use(N_x/2+1:N_x) = force_next_corr(N_x/2+1:N_x);
            end

            
            % Update velocity left
            if method_left <= 2 || method_left == 5
                % (unused) curr force not corrected in second order
                % next force corrected in first order
                v(1:N_x/2,n+1) = update_velocity_FDTD(data_left, p(1:N_x/2,n+1), p(1:N_x/2,n), p(1:N_x/2,n-1), force_use(1:N_x/2), v(1:N_x/2,n), override_order);
            elseif method_left >= 3
                % curr force not corrected in second order
                % next force corrected in first order
                v(1:N_x/2,n+1) = update_velocity_Fourier(data_left, p(1:N_x/2,n+1), p(1:N_x/2,n), p(1:N_x/2,n-1), force_use(1:N_x/2), v(1:N_x/2,n), override_order);
            end
            
            % Update velocity right
            if method_right <= 2 || method_right == 5
                v(N_x/2+1:N_x,n+1) = update_velocity_FDTD(data_right, p(N_x/2+1:N_x,n+1), p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_use(N_x/2+1:N_x), v(N_x/2+1:N_x,n), override_order);
            elseif method_right >= 3
                v(N_x/2+1:N_x,n+1) = update_velocity_Fourier(data_right, p(N_x/2+1:N_x,n+1), p(N_x/2+1:N_x,n), p(N_x/2+1:N_x,n-1), force_use(N_x/2+1:N_x), v(N_x/2+1:N_x,n), override_order);
            end
            
            % Post-merge first order
            % using r^{n+1}, on v^{n+1}
            if merge_left == 2
                if temp_order_left == 1 || override_order == true
                    v(1:N_x/2,n+1) = v(1:N_x/2,n+1) + transmittivity^2 * dt * residual(1:N_x/2) / (1 + 2*dt*alpha_abs);
                end
            end

            if merge_right == 2
                if temp_order_right == 1 || override_order == true
                    v(N_x/2+1:N_x,n+1) = v(N_x/2+1:N_x,n+1) + transmittivity^2 * dt * residual(N_x/2+1:N_x) / (1 + 2*dt*alpha_abs);
                end
            end

        else % no DD

            % Update force
            force_now = force_next;

            % Update pressure
            if method_left <= 2 || method_left == 5
                p(:,n+1) = update_pressure_FDTD(data_left, p(:,n), p(:,n-1), force_now(:), v(:,n), override_order);
            elseif method_left >= 3
                p(:,n+1) = update_pressure_Fourier(data_left, p(:,n), p(:,n-1), force_now(:), v(:,n), override_order);
            end

            % Impose next force
            force_next = force(:,n+1);

            % Artificial dissipation for stability
            p(:,n+1) = diss * p(:,n+1);

            % Update velocity
            if method_left <= 2 || method_left == 5
                v(:,n+1) = update_velocity_FDTD(data_left, p(:,n+1), p(:,n), p(:,n-1), force_next(:), v(:,n), override_order);
            elseif method_left >= 3
                v(:,n+1) = update_velocity_Fourier(data_left, p(:,n+1), p(:,n), p(:,n-1), force_next(:), v(:,n), override_order);
            end

        end
        
        info_str = ['Instant [s]: ' num2str(n*dt, '%4.3f') ' / ' ...
                num2str(len_t, '%4.3f') ' ( ' num2str(n/N_t*100, '%4.1f') '% )'];

        if debug == true
            if mod(n-1,10) == 1
	            plot_snapshot(x_axis,len_x,p(:,n+1),v(:,n+1),f,db_plot);
            end
            figure(f);
            sgtitle(info_str);
        else
            clc;
            disp(info_str);
        end
    
    end

    % Discarding stub time instant
    p = p(:,2:end);
    v = v(:,2:end);
    
end

