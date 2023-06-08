function simulation_parameters = get_simulation_parameters()
    msg4 = "Use domain decomposition?";
    opts4 = ["Yes" "No"];
    choice4 = menu(msg4, opts4);

    DD = choice4 == 1;

    choice_merge_left = 0; % placeholder
    choice_merge_right = 0; % placeholder

    if DD
        msg_merge_left = "Choose the merge approach for left";
        opts_merge = ["Pre-merge" "Post-merge"];
        choice_merge_left = menu(msg_merge_left, opts_merge);

        msg_merge_right = "Choose the merge approach for right";
        choice_merge_right = menu(msg_merge_right, opts_merge);

        msg_method_left = "Choose the update method for left";
        opts_method = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
        choice_method_left = menu(msg_method_left, opts_method);
    
        msg_method_right = "Choose the update method for right";
        choice_method_right = menu(msg_method_right, opts_method);
    
        if (choice_method_left == 2) || (choice_method_left == 4)
            order_left = 1;
        else
            order_left = 2;
        end
    
        if (choice_method_right == 2) || (choice_method_right == 4)
            order_right = 1;
        else
            order_right = 2;
        end
    else
        msg_method_left = "Choose the update method";
        opts_method = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
        choice_method_left = menu(msg_method_left, opts_method);
        choice_method_right = choice_method_left; % placeholder
    
        if (choice_method_left == 2) || (choice_method_left == 4)
            order_left = 1;
            order_right = order_left;
        else
            order_left = 2;
            order_right = order_left;
        end
    end

    simulation_parameters.merge_left = choice_merge_left;
    simulation_parameters.merge_right = choice_merge_right;
    simulation_parameters.method_left = choice_method_left;
    simulation_parameters.method_right = choice_method_right;
    simulation_parameters.order_left = order_left;
    simulation_parameters.order_right = order_right;
    simulation_parameters.DD = DD;
end
