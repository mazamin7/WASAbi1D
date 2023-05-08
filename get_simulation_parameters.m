function simulation_parameters = get_simulation_parameters()
    msg4 = "Use domain decomposition?";
    opts4 = ["Yes" "No"];
    choice4 = menu(msg4, opts4);

    DD = choice4 == 1;

    choice = 0; % placeholder

    if DD
        msg = "Choose the merge approach";
        opts = ["Pre 7 points" "Post 7 points"];
        choice = menu(msg, opts);

        msg2 = "Choose the update method for left";
        opts2 = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
        choice2 = menu(msg2, opts2);
    
        msg3 = "Choose the update method for right";
        choice3 = menu(msg3, opts2);
    
        if (choice2 == 2) || (choice2 == 4)
            order_left = 1;
        else
            order_left = 2;
        end
    
        if (choice3 == 2) || (choice3 == 4)
            order_right = 1;
        else
            order_right = 2;
        end
    
        assert(order_left == order_right, 'The order of left and right method must be the same');
        % otherwise domain decomposition doesn't work
    
        order = order_left;
    else
        msg2 = "Choose the update method";
        opts2 = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
        choice2 = menu(msg2, opts2);
        choice3 = choice2; % placeholder
    
        if (choice2 == 2) || (choice2 == 4)
            order = 1;
        else
            order = 2;
        end
    end

    simulation_parameters.merge = choice;
    simulation_parameters.method_left = choice2;
    simulation_parameters.method_right = choice3;
    simulation_parameters.order = order;
    simulation_parameters.DD = DD;
end
