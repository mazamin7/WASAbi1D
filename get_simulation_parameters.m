function simulation_parameters = get_simulation_parameters()
    msg = "Choose the merge approach";
    opts = ["Pre 7 points" "Post 7 points"];
    choice = menu(msg, opts);

    msg2 = "Choose the update for left";
    opts2 = ["FDTD 2ord" "FDTD 1ord" "Fourier 2ord" "Fourier 1ord" "PML"];
    choice2 = menu(msg2, opts2);

    msg3 = "Choose the update for right";
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

    simulation_parameters.merge = choice;
    simulation_parameters.method_left = choice2;
    simulation_parameters.method_right = choice3;
    simulation_parameters.order = order;
end
