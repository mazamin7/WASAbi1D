function save_plots(test_case_data, simulation_parameters, dt, dh, fig_p, fig_v)

    % Extracting test case data
    test_case = test_case_data.test_case;
    alpha_abs = test_case_data.alpha_abs;
    bc_left = test_case_data.bc_left;
    bc_right = test_case_data.bc_right;

    % Extracting simulation parameters
    merge = simulation_parameters.merge;
    method_left = simulation_parameters.method_left;
    method_right = simulation_parameters.method_right;
    DD = simulation_parameters.DD;
    
    damped = alpha_abs ~= 0;

    % Save figures as images
    if merge == 1
        merge_str = 'pre_merge';
    elseif merge == 2
        merge_str = 'post_merge';
    elseif DD == false
        merge_str = 'noDD';
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
    
    if damped == true
        damping_str = 'damped';
    else
        damping_str = 'undamped';
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
    foldername = sprintf('test=%s__dh=%.4f_dt=%.4f__%s_%s__left=%s_%s__right=%s_%s', num2str(test_case), dh, dt, merge_str, damping_str, method_left_str, bc_left_str, method_right_str, bc_right_str);

    % Create images folder if not exists
    if ~exist('images', 'dir')
        mkdir('images')
    end
    
    mkdir(fullfile('images', foldername));
    
    % Save figures as images in the folder
    filename_p_top = fullfile('images', foldername, 'pressure_top.png');
    filename_p_3d_a = fullfile('images', foldername, 'pressure_3d_a.png');
    filename_p_3d_b = fullfile('images', foldername, 'pressure_3d_b.png');
    filename_v_top = fullfile('images', foldername, 'velocity_top.png');
    filename_v_3d_a = fullfile('images', foldername, 'velocity_3d_a.png');
    filename_v_3d_b = fullfile('images', foldername, 'velocity_3d_b.png');
        
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
        % 3D view A for pressure
        figure(fig_p)
        view(3)
        saveas(fig_p, filename_p_3d_a);
        
        % 3D view A for velocity
        figure(fig_v)
        view(3)
        saveas(gcf, filename_v_3d_a);


        % 3D view B for pressure
        figure(fig_p)
        view(30,30)
        saveas(fig_p, filename_p_3d_b);
        
        % 3D view B for velocity
        figure(fig_v)
        view(30,30)
        saveas(gcf, filename_v_3d_b);
    end

end