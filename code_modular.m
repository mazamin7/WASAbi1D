clear all, close all, clc;

% User menu
msg = "Choose the merge approach";
opts = ["Pre 7 points" "Post 7 points"];
choice = menu(msg, opts);

msg2 = "Choose the update for left";
opts2 = ["FDTD" "Fourier" "FEM" "PML"];
choice2 = menu(msg2, opts2);

msg3 = "Choose the update for right";
choice3 = menu(msg3, opts2);

if choice2 < 4
    msg4 = "Is left damped?";
    opts4 = ["Yes" "No"];
    choice4 = menu(msg4, opts4);
else
    choice4 = 1;
end

if choice3 < 4
    msg5 = "Is right damped?";
    choice5 = menu(msg5, opts4);
else
    choice5 = 1;
end

opts6 = ["N" "D"];

if choice2 ~= 4
    msg6 = "Choose boundary condition for left";
    choice6 = menu(msg6, opts6);

    boundCondLeft = opts6(choice6);
else
    boundCondLeft = "N";
end

if choice3 ~= 4
    msg7 = "Choose boundary condition for right";
    choice7 = menu(msg7, opts6);

    boundCondRight = opts6(choice7);
else
    boundCondRight = "N";
end

% Decide whether FDTD treats boundaries explicitly or not
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

q_next_dct_left = zeros(N/2,1);
q_next_dct_right = zeros(N/2,1);

% Imposing initial conditions
pulse_width = 1/2^4;
pulse_pos = 1/2;

pulse_width_x = floor(pulse_width * N);
pulse_pos_x = floor(pulse_pos * N);

pulse_axis = 1:pulse_width_x;
pulse = 1/2 - 1/2 * cos(2*pi*pulse_axis/pulse_width_x);

% p_curr(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;
% p_prev(pulse_pos_x-pulse_width_x/2+1:pulse_pos_x+pulse_width_x/2) = pulse;

% Building pre/post-merge matrix
alpha = 1/90;
beta = -3/20;
gamma = 3/2;
delta = -49/18;

C = sparse(N,N);

if shiftLeft == false
    C(N/2-2,N/2:N/2+1) = [-alpha, alpha];
    C(N/2-1,N/2-1:N/2+2) = [-alpha, -beta, beta, alpha];
    C(N/2,N/2-2:N/2+3) = [-alpha, -beta, -gamma, gamma, beta, alpha];
else
    C(N/2-2,N/2-1:N/2+2) = [-alpha, 0, 0, alpha];
    C(N/2-1,N/2-2:N/2+3) = [-alpha, -beta, 0, 0, beta, alpha];
    C(N/2,N/2-3:N/2+4) = [-alpha, -beta, -gamma, 0, 0, gamma, beta, alpha];
end

if shiftRight == false
    C(N/2+1,N/2-2:N/2+3) = -[-alpha, -beta, -gamma, gamma, beta, alpha];
    C(N/2+2,N/2-1:N/2+2) = -[-alpha, -beta, beta, alpha];
    C(N/2+3,N/2:N/2+1) = -[-alpha, alpha];
else
    C(N/2+1,N/2-3:N/2+4) = -[-alpha, -beta, -gamma, 0, 0, gamma, beta, alpha];
    C(N/2+2,N/2-2:N/2+3) = -[-alpha, -beta, 0, 0, beta, alpha];
    C(N/2+3,N/2-1:N/2+2) = -[-alpha, 0, 0, alpha];
end

C_leftBC = sparse(N,N);
C_rightBC = sparse(N,N);

C_leftBC(1,1:3) = -[gamma beta alpha];
C_leftBC(2,1:2) = -[beta alpha];
C_leftBC(3,1) = -alpha;

C_rightBC(N-2,N) = -alpha;
C_rightBC(N-1,N-1:N) = -[beta alpha];
C_rightBC(N,N-2:N) = -[gamma beta alpha];

% Initializing update methods
if choice2 == 1 || choice2 == 4
    FDTD_data_left = init_FDTD(N/2, c, dt, dh, choice4 == 1, alpha_abs, explicitBoundariesFDTD == true, boundCondLeft, "N", choice2 == 4);
elseif choice2 == 2
    Fourier_data_left = init_Fourier(N/2, c, dt, dh, choice4 == 1, alpha_abs);
else
    FEM_data_left = init_FEM(N/2, c, dt, dh, choice4 == 1, alpha_abs, boundCondLeft, "N");
end

if choice3 == 1 || choice3 == 4
    FDTD_data_right = init_FDTD(N/2, c, dt, dh, choice5 == 1, alpha_abs, explicitBoundariesFDTD == true, "N", boundCondRight, choice3 == 4);
elseif choice3 == 2
    Fourier_data_right = init_Fourier(N/2, c, dt, dh, choice5 == 1, alpha_abs);
else
    FEM_data_right = init_FEM(N/2, c, dt, dh, choice4 == 1, alpha_abs, "N", boundCondRight);
end

% Simulation loop
for n = 1:dur_samples

    force = zeros(N,1);
    % force(floor(N/2)) = 1000*sin(2*pi*6.2832*n*dt);

    g1_dirichlet = 0.2*sin(2*pi*4*n*dt);
    g2_dirichlet = g1_dirichlet;

    % Pre-merge
    if choice == 1
        force = force + (c / dh)^2 * C * p_curr;
    end
    
    % Update left
    if choice2 == 1 || choice2 == 4
        p_next(1:N/2) = update_FDTD(FDTD_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2));
    elseif choice2 == 2
        [p_next(1:N/2),q_next_dct_left] = update_Fourier(Fourier_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2), q_next_dct_left);
    else
        p_next(1:N/2) = update_FEM(FEM_data_left, p_curr(1:N/2), p_prev(1:N/2), force(1:N/2));
    end
    
    % Update right
    if choice3 == 1 || choice3 == 4
        p_next(N/2+1:N) = update_FDTD(FDTD_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N));
    elseif choice3 == 2
        [p_next(N/2+1:N),q_next_dct_right] = update_Fourier(Fourier_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N), q_next_dct_right);
    else
        p_next(N/2+1:N) = update_FEM(FEM_data_right, p_curr(N/2+1:N), p_prev(N/2+1:N), force(N/2+1:N));
    end

    % B.C. for Fourier
    if choice2 == 2
        if choice6 == 2

            % OLD METHOD
            % p_next = p_next + (c * dt / dh)^2 * C_leftBC * p_curr;

            % NEW METHOD
            p_next(2) = g1_dirichlet;

        end
    end

    if choice3 == 2
        if choice7 == 2

            % OLD METHOD
            % p_next = p_next + (c * dt / dh)^2 * C_rightBC * p_curr;

            % NEW METHOD
            p_next(N-1) = g2_dirichlet;

        end
    end
    
    % Post-merge
    if choice == 2
        p_next = p_next + (c * dt / dh)^2 * C * p_curr;
    end

    % Fixing interfaces
    if shiftLeft == true && shiftRight == true
        p_next(N/2+1) = 0.5*(p_next(N/2) + p_next(N/2+1));
        p_next(N/2) = p_next(N/2+1);
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

    pause(0.1);

end