function save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v)

    % Extracting test case data
    test_case = test_case_data.test_case;
    alpha_abs_left = test_case_data.alpha_abs_left;
    alpha_abs_right = test_case_data.alpha_abs_right;
    bc_left = test_case_data.bc_left;
    bc_right = test_case_data.bc_right;

    % Extracting simulation parameters
    merge = simulation_parameters.merge;
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    
    if alpha_abs_left == 0
        left_damped = false;
    else
        left_damped = true;
    end
    
    if alpha_abs_right == 0
        right_damped = false;
    else
        right_damped = true;
    end

    % Save figures as images
    if merge == 1
        merge_str = 'pre_merge';
    elseif merge == 2
        merge_str = 'post_merge';
    end
    
    if method_left == 1
        method_left_str = 'FDTD_2ord';
    elseif method_left == 2
        method_left_str = 'FDTD_1ord';
    elseif method_left == 3
        method_left_str = 'Fourier_2ord';
    elseif method_left == 4
        method_left_str = 'Fourier_1ord';
    elseif method_left == 5
        method_left_str = 'PML';
    end
    
    if method_right == 1
        method_right_str = 'FDTD_2ord';
    elseif method_right == 2
        method_right_str = 'FDTD_1ord';
    elseif method_right == 3
        method_right_str = 'Fourier_2ord';
    elseif method_right == 4
        method_right_str = 'Fourier_1ord';
    elseif method_right == 5
        method_right_str = 'PML';
    end
    
    if left_damped == true
        damping_left_str = 'damped';
    else
        damping_left_str = 'undamped';
    end
    
    if right_damped == true
        damping_right_str = 'damped';
    else
        damping_right_str = 'undamped';
    end
    
    if bc_left == 'D'
        bc_left_str = 'Dirichlet';
    elseif bc_left == 'N'
        bc_left_str = 'Neumann';
    end
    
    if bc_right == 'D'
        bc_right_str = 'Dirichlet';
    elseif bc_right == 'N'
        bc_right_str = 'Neumann';
    end
    
    % Create folder with current filename
    foldername = sprintf('test=%s__dh=%.2f_dt=%.2f__%s__left=%s_%s_%s__right=%s_%s_%s', num2str(test_case), dh, dt, merge_str, method_left_str, damping_left_str, bc_left_str, method_right_str, damping_right_str, bc_right_str);

    % Create images folder if not exists
    if ~exist('images', 'dir')
        mkdir('images')
    end
    
    mkdir(fullfile('images', foldername));
    
    % Save figures as images in the folder
    filename_p_top = fullfile('images', foldername, 'pressure_top.png');
    filename_p_3d = fullfile('images', foldername, 'pressure_3d.png');
    filename_v_top = fullfile('images', foldername, 'velocity_top.png');
    filename_v_3d = fullfile('images', foldername, 'velocity_3d.png');
        
    view_top = true;
    view_3d = true;
    
    if view_top
        % Top view for pressure
        figure(fig_p)
        view(0,90)
        saveas(gcf, filename_p_top);
        
        % Top view for velocity
        figure(fig_v)
        view(0,90)
        saveas(gcf, filename_v_top);
    end
    
    if view_3d
        % 3D view for pressure
        figure(fig_p)
        view(3)
        saveas(fig_p, filename_p_3d);
        
        % 3D view for velocity
        figure(fig_v)
        view(3)
        saveas(gcf, filename_v_3d);
    end

end