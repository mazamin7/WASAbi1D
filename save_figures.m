function save_plots(choice, choice2, choice3, left_damped, right_damped, bc_left, bc_right, dh, dt)

    % Save figures as images
    if choice == 1
        merge_str = 'pre_merge';
    elseif choice == 2
        merge_str = 'post_merge';
    end
    
    if choice2 == 1
        method_left_str = 'FDTD_2ord';
    elseif choice2 == 2
        method_left_str = 'FDTD_1ord';
    elseif choice2 == 3
        method_left_str = 'Fourier_2ord';
    elseif choice2 == 4
        method_left_str = 'Fourier_1ord';
    elseif choice2 == 5
        method_left_str = 'PML';
    end
    
    if choice3 == 1
        method_right_str = 'FDTD_2ord';
    elseif choice3 == 2
        method_right_str = 'FDTD_1ord';
    elseif choice3 == 3
        method_right_str = 'Fourier_2ord';
    elseif choice3 == 4
        method_right_str = 'Fourier_1ord';
    elseif choice3 == 5
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
    foldername = sprintf('dh=%.2f_dt=%.2f__%s__left=%s_%s_%s__right=%s_%s_%s', dh, dt, merge_str, method_left_str, damping_left_str, bc_left_str, method_right_str, damping_right_str, bc_right_str);
    
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
        figure(3)
        view(0,90)
        saveas(gcf, filename_p_top);
        
        % Top view for velocity
        figure(4)
        view(0,90)
        saveas(gcf, filename_v_top);
    end
    
    if view_3d
        % 3D view for pressure
        figure(3)
        view(3)
        saveas(gcf, filename_p_3d);
        
        % 3D view for velocity
        figure(4)
        view(3)
        saveas(gcf, filename_v_3d);
    end

end