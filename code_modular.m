clear all, close all, clc;

% User menu
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

% Simulation parameters
dh = 0.1;
dt = 0.002;

[len_x, len_t, c0, alpha_abs_left, alpha_abs_right, transmittivity, ...
    bc_left, bc_right, force_time_fun, force_envelope, ...
    g1_time_fun, g2_time_fun] = get_test_case(1);

% Checking compatibility between simulation parameters and test case
assert(dt < dh / 2 / c0);

assert(~((choice2 == 3 || choice2 == 4) && (bc_left == "D" || bc_right == "D")), ...
    'Boundary conditions must be NEUMANN homogeneous with Fourier method')

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

assert(~(choice2 == 3 && left_damped));
assert(~(choice3 == 3 && right_damped));



% Knowing simulation pars and test case, initialize simulation variables

% Defining time and space axes
N_t = floor(len_t / dt);
N_x = floor(len_x / dh);

x_axis = linspace(0,len_x,N_x);
t_axis = linspace(0,len_t,N_t);

% Building residue matrix
C = get_residue_matrix(N_x, 6);

% Initialize solution, force, and boundary conditions data
p = zeros(N_x,N_t);
v = zeros(N_x,N_t);
force = zeros(N_x,N_t);
g1 = zeros(N_t,1);
g2 = zeros(N_t,1);

for n = 1:N_t
    force(:,n) = force_envelope(x_axis) * force_time_fun(t_axis(n));
    g1(n) = g1_time_fun(t_axis(n));
    g2(n) = g2_time_fun(t_axis(n));
end

if (choice2 == 3 || choice2 == 4 || choice3 == 3 || choice3 == 4)
    assert(all(g1 == 0), 'Boundary conditions must be Neumann HOMOGENEOUS with Fourier method');
    assert(all(g2 == 0), 'Boundary conditions must be Neumann HOMOGENEOUS with Fourier method');
end

% Initializing update methods
if choice2 <= 2 || choice2 == 5
    data_left = init_FDTD(len_x/2, c0, dt, dh, alpha_abs_left, bc_left, "N", choice2 == 5, order_left);
elseif choice2 >= 3
    data_left = init_Fourier(len_x/2, c0, dt, dh, order_left, alpha_abs_left);
end

if choice3 <= 2 || choice3 == 5
    data_right = init_FDTD(len_x/2, c0, dt, dh, alpha_abs_right, "N", bc_right, choice3 == 5, order_right);
elseif choice3 >= 3
    data_right = init_Fourier(len_x/2, c0, dt, dh, order_right, alpha_abs_right);
end

% Simulation loop
for n = 3:N_t

    % Residual calculation
    residual = (c0 / dh)^2 * C * p(:,n-1);

    % Pre-merge
    if choice == 1
        force_now = force(:,n-1) + transmittivity^2 * residual;
    else
        force_now = force(:,n-1);
    end
    
    % Update left
    if choice2 <= 2 || choice2 == 5
        [p(1:N_x/2,n),v(1:N_x/2,n)] = update_FDTD(data_left, p(1:N_x/2,n-1), p(1:N_x/2,n-2), force_now(1:N_x/2), v(1:N_x/2,n-1), v(1:N_x/2,n-2), g1(n-1), 0);
    elseif choice2 >= 3
        [p(1:N_x/2,n),v(1:N_x/2,n)] = update_Fourier(data_left, p(1:N_x/2,n-1), p(1:N_x/2,n-2), force_now(1:N_x/2), v(1:N_x/2,n-1), v(1:N_x/2,n-2));
    end
    
    % Update right
    if choice3 <= 2 || choice3 == 5
        [p(N_x/2+1:N_x,n),v(N_x/2+1:N_x,n)] = update_FDTD(data_right, p(N_x/2+1:N_x,n-1), p(N_x/2+1:N_x,n-2), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n-1), v(N_x/2+1:N_x,n-2), 0, g2(n-1));
    elseif choice3 >= 3
        [p(N_x/2+1:N_x,n),v(N_x/2+1:N_x,n)] = update_Fourier(data_right, p(N_x/2+1:N_x,n-1), p(N_x/2+1:N_x,n-2), force_now(N_x/2+1:N_x), v(N_x/2+1:N_x,n-1), v(N_x/2+1:N_x,n-2));
    end

    % Post-merge
    if choice == 2
        if order_left == 1
            v(:,n) = v(:,n) + transmittivity^2 * 2*dt * residual;
        elseif order_left == 2
            p(:,n) = p(:,n) + transmittivity^2 * dt*dt * residual;
        end
    end
    
    % Plot
    f = figure(2);
    f.Position = [100, 100, 1200, 700];
    sgtitle(['instant [s]: ' num2str(n*dt, '%4.3f') ' / ' ...
        num2str(len_t, '%4.3f') ' ( ' num2str(n/N_t*100, '%4.1f') '% )']);

    % Plot p
    subplot(2,1,1);
    plot(x_axis, p(:,n));
    title('Pressure');
    xlim([0,len_x]);
    ylim([-1,1]*2e-1);

    % Plot v
    subplot(2,1,2);
    plot(x_axis, v(:,n));
    title('Velocity');
    xlim([0,len_x]);
    ylim([-c0,c0]*5e-1);

end

% Plot p
figure(3);
surf(t_axis, x_axis, p,'EdgeColor','none','FaceColor','interp');
xlabel('Time [s]');
ylabel('Space [m]');
zlabel('Pressure');
title('Pressure Solution');
view(0, 90);  % set view to show from the top

% Plot v
figure(4);
surf(t_axis, x_axis, v,'EdgeColor','none','FaceColor','interp');
xlabel('Time [s]');
ylabel('Space [m]');
zlabel('Velocity');
title('Velocity Solution');
view(0, 90);  % set view to show from the top

% Save simulation as figures and animation
save_plots(choice, choice2, choice3, left_damped, right_damped, bc_left, bc_right, dh, dt);







