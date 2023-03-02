clear all, close all, clc;

% User menu
msg = "Choose the merge approach";
opts = ["Pre 7 points" "Post 7 points" "Borrel 3 points" "Borrel 7 points"];
choice = menu(msg, opts);

msg2 = "Choose the update for left";
opts2 = ["FDTD" "Fourier" "FEM"];
choice2 = menu(msg2, opts2);

msg3 = "Choose the update for right";
choice3 = menu(msg3, opts2);

msg4 = "Is left damped?";
opts4 = ["Yes" "No"];
choice4 = menu(msg4, opts4);

msg5 = "Is right damped?";
choice5 = menu(msg5, opts4);

opts6 = ["N" "D"];

if choice2 ~= 2
    msg6 = "Choose boundary condition for left";
    choice6 = menu(msg6, opts6);

    boundCondLeft = opts6(choice6);
end

if choice3 ~= 2
    msg7 = "Choose boundary condition for right";
    choice7 = menu(msg7, opts6);

    boundCondRight = opts6(choice7);
end

% Decide whether FDTD treats boundaries explicitly or not
% Only configurable with Borrel-merge
explicitBoundariesFDTD = false; % WORKS BETTER IF SET TO FALSE
% false -> boundary at N/2
% true -> boundary at N

% Shift values for Borrel-merge (treat explicitly boundary)
shiftLeft = choice2 == 3 || (choice2 == 1 && explicitBoundariesFDTD == true);
shiftRight = choice3 == 3 || (choice3 == 1 && explicitBoundariesFDTD == true);

% Simulation parameters
N = 2^7;
dt = 1/(2*N);
c = 1;

alpha_abs = 10; % Absorption coefficient

len = 1; % Domain length
dur = 50; % Simulation duration

% Defining time and space axis
dur_samples = floor(dur / dt);
dh = len/(N-1);
x_axis = 1:N;

% Checking parameters validity
assert(dt <= dh / sqrt(3) / c)

% Initializing solution data
p_prev = zeros(N,1);
p_curr = p_prev * 0;
p_next = p_prev * 0;

% Imposing initial conditions
pulse_width = 1/2^3;
pulse_pos = 3/4;

pulse_width_x = floor(pulse_width * N);
pulse_pos_x = floor(pulse_pos * N);

pulse_axis = 1:pulse_width_x;
pulse = 1/2 - 1/2 * cos(2*pi*pulse_axis/pulse_width_x);

p_curr(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;
p_prev(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;

% Building pre/post-merge matrix
alpha = 1/90;
beta = -3/20;
gamma = 3/2;
delta = -49/18;

C = sparse(N,N);

C(N/2-2,N/2:N/2+1) = [-alpha, alpha];
C(N/2-1,N/2-1:N/2+2) = [-alpha, -beta, beta, alpha];
C(N/2,N/2-2:N/2+3) = [-alpha, -beta, -gamma, gamma, beta, alpha];
C(N/2+1,N/2-2:N/2+3) = -C(N/2,N/2-2:N/2+3);
C(N/2+2,N/2-1:N/2+2) = -C(N/2-1,N/2-1:N/2+2);
C(N/2+3,N/2:N/2+1) = -C(N/2-2,N/2:N/2+1);

% Initializing update methods
if choice2 == 1
    FDTD_data_left = init_FDTD(N/2, c, dt, dh, choice4 == 1, alpha_abs, choice > 2 && explicitBoundariesFDTD == true, boundCondLeft, "N");
elseif choice2 == 2
    Fourier_data_left = init_Fourier(N/2, c, dt, dh, choice4 == 1, alpha_abs);
else
    FEM_data_left = init_FEM(N/2, c, dt, dh, choice4 == 1, alpha_abs, boundCondLeft, "N");
end

if choice3 == 1
    FDTD_data_right = init_FDTD(N/2, c, dt, dh, choice5 == 1, alpha_abs, choice > 2 && explicitBoundariesFDTD == true, "N", boundCondRight);
elseif choice3 == 2
    Fourier_data_right = init_Fourier(N/2, c, dt, dh, choice5 == 1, alpha_abs);
else
    FEM_data_right = init_FEM(N/2, c, dt, dh, choice4 == 1, alpha_abs, "N", boundCondRight);
end

% Simulation loop
for n = 1:dur_samples

    force = zeros(N,1);
    % force(floor(N/2)) = 1000*sin(2*pi*6.2832*n*dt);

    % Pre-merge
    if choice == 1
        force = force + (c / dh)^2 * C * p_curr;
    end
    
    % Update left
    if choice2 == 1
        p_next(1:N/2) = update_FDTD(FDTD_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2));
    elseif choice2 == 2
        p_next(1:N/2) = update_Fourier(Fourier_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2));
    else
        p_next(1:N/2) = update_FEM(FEM_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2));
    end
    
    % Update right
    if choice3 == 1
        p_next(N/2+1:N) = update_FDTD(FDTD_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N));
    elseif choice3 == 2
        p_next(N/2+1:N) = update_Fourier(Fourier_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N));
    else
        p_next(N/2+1:N) = update_FEM(FEM_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N));
    end
    
    % Post-merge
    if choice == 2
        p_next = p_next + (c * dt / dh)^2 * C * p_curr;
    elseif choice == 3
        % Compute residual parts
        if shiftLeft == false
            res1 = c^2 * dt^2 / dh^2 * p_curr(N/2);
        else
            res1 = c^2 * dt^2 / dh^2 * p_curr(N/2-1);
        end

        if shiftRight == false
            res2 = c^2 * dt^2 / dh^2 * p_curr(N/2+1);
        else
            res2 = c^2 * dt^2 / dh^2 * p_curr(N/2+2);
        end

        % Remove the residual part and transfer the removed residual part to the other domain
        p_next(N/2) = p_next(N/2) - res1 + res2;
        p_next(N/2+1) = p_next(N/2+1) - res2 + res1;
    elseif choice == 4
        % Compute residual parts
        if shiftLeft == false
            res1 = c^2 * dt^2 / dh^2 * [alpha, beta, gamma] * p_curr(N/2-2:N/2);
        else
            res1 = c^2 * dt^2 / dh^2 * [alpha, beta, gamma] * p_curr(N/2-3:N/2-1);
        end

        if shiftRight == false
            res2 = c^2 * dt^2 / dh^2 * [gamma, beta, alpha] * p_curr(N/2+1:N/2+3);
        else
            res2 = c^2 * dt^2 / dh^2 * [gamma, beta, alpha] * p_curr(N/2+2:N/2+4);
        end

        % Remove the residual part and transfer the removed residual part to the other domain
        p_next(N/2) = p_next(N/2) - res1 + res2;
        p_next(N/2+1) = p_next(N/2+1) - res2 + res1;
    end
    
    % Update
    p_prev = p_curr;
    p_curr = p_next;
    
    % Plot
    f = figure(1);
    f.Position = [100, 100, 1500, 400];
    plot(x_axis*dh, full(p_next));
    xlim([0,len])
    ylim([-1,1]);

    sgtitle(['instant [s]: ' num2str((n+1)*dt, '%4.3f') ' / ' num2str(dur, '%4.3f') ' ( ' num2str((n+1)/dur_samples*100, '%4.1f') '% )']);

end