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

assert(choice2 ~= 3 && choice3 ~= 3, 'FEM not supported');

% Simulation parameters
c0 = 343;
len_x = 20; % Domain length
T_sec = 0.2; % Simulation duration
alpha_abs_left = 0; % Absorption coefficient
alpha_abs_right = 50;
exact_damping = true; % Use exact damping
dh = 1/2^5;
transmittivity = 1; % Transmittance of the middle boundary

dt = dh / 2 / c0;

% Boundary conditions
bc_left = "N";
bc_right = "N";

% Is partition damped?
left_damped = true;
right_damped = true;

assert(exact_damping == false || left_damped == right_damped);

% Source
freq_source = 300;
source_fun = @(n) 8e5 * sin(2*pi*freq_source*n*dt) * (n*dt <= 1/freq_source);
source_pos_ratio_x = 2/10;
sigma = len_x/40;       % standard deviation of force spatial envelope (Gaussian)
mach_x = 0;

% Checking parameters validity
assert(dt <= dh / sqrt(3) / c0)

% Defining time and space axis
N_t = floor(T_sec / dt);
N_x = floor(len_x / dh);

% Building residue matrix
C = get_residue_matrix(N_x, 8);

% Initialize solution data
p_prev = zeros(N_x,1);
p_curr = zeros(N_x,1);
p_next = zeros(N_x,1);
q_next = zeros(N_x,1); % for exact viscous damping
q_curr = zeros(N_x,1); % for exact viscous damping

% Defining force spatial envelope
x_axis = linspace(0,len_x,N_x);
force_pos = len_x * source_pos_ratio_x;
force_envelope = @(x,mu) 1/(sigma * sqrt(2 * pi)) * exp(-(x-mu).^2/(2*sigma^2)); % Gaussian function

% Initializing update methods
if choice2 == 1 || choice2 == 4
    FDTD_data_left = init_FDTD(len_x/2, c0, dt, dh, left_damped, alpha_abs, bc_left, "N", choice2 == 4);
elseif choice2 == 2
    Fourier_data_left = init_Fourier(len_x/2, c0, dt, dh, left_damped, alpha_abs_left);
else
    FEM_data_left = init_FEM(len_x/2, c0, dt, dh, left_damped, alpha_abs, bc_left, "N");
end

if choice3 == 1 || choice3 == 4
    FDTD_data_right = init_FDTD(len_x/2, c0, dt, dh, right_damped, alpha_abs, "N", bc_right, choice3 == 4);
elseif choice3 == 2
    Fourier_data_right = init_Fourier(len_x/2, c0, dt, dh, right_damped, alpha_abs_right);
else
    FEM_data_right = init_FEM(len_x/2, c0, dt, dh, right_damped, alpha_abs, "N", bc_right);
end

% Simulation loop
for n = 1:N_t

    % Update force position
    force_pos = force_pos + mach_x * c0 * dt;

    % Impose force
    force_envelope_temp = force_envelope(x_axis,force_pos);
    force = source_fun(n) * force_envelope_temp';

%     % Plot force envelope
%     figure(1);
%     plot(x_axis, force_envelope_temp)
%     xlabel('Position')
%     ylabel('Intensity')
%     title('Gaussian Force Envelope')

    % Impose non-homogeneous b.c.
    g1 = 0; % 1/2*0.3*sin(2*pi*4*n*dt) * (n <= 1/4 / dt);
    g2 = 0; % g1;
    % FOR NOW, ONLY FDTD SUPPORTS NON-HOMOGENEOUS DIRICHLET/NEUMANN

    % Pre-merge
    if choice == 1
        force = force + transmittivity^2 * (c0 / dh)^2 * C * p_curr;
    end
    
    % Update left
    if choice2 == 1 || choice2 == 4
        p_next(1:N_x/2) = update_FDTD(FDTD_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2), g1, 0);
    elseif choice2 == 2
        [p_next(1:N_x/2),q_next(1:N_x/2)] = update_Fourier(Fourier_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2), q_curr(1:N_x/2), exact_damping);
    else
        p_next(1:N_x/2) = update_FEM(FEM_data_left, p_curr(1:N_x/2), p_prev(1:N_x/2), force(1:N_x/2));
    end
    
    % Update right
    if choice3 == 1 || choice3 == 4
        p_next(N_x/2+1:N_x) = update_FDTD(FDTD_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x), 0, g2);
    elseif choice3 == 2
        [p_next(N_x/2+1:N_x),q_next(N_x/2+1:N_x)] = update_Fourier(Fourier_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x), q_curr(N_x/2+1:N_x), exact_damping);
    else
        p_next(N_x/2+1:N_x) = update_FEM(FEM_data_right, p_curr(N_x/2+1:N_x), p_prev(N_x/2+1:N_x), force(N_x/2+1:N_x));
    end

    % Post-merge
    if choice == 2
        p_next = p_next + transmittivity * (c0 * dt / dh)^2 * C * p_curr;
        q_next = q_next + (c0 * dt / dh)^2 * C * q_curr;
    end

    % Update
    p_prev = p_curr;
    p_curr = p_next;
    q_curr = q_next;
    
    % Plot
    f = figure(2);
    f.Position = [100, 100, 1500, 400];
    plot(x_axis, p_next);
    xlim([0,len_x]);
    ylim([-1,1]);

    sgtitle(['instant [s]: ' num2str((n+1)*dt, '%4.3f') ' / ' num2str(T_sec, '%4.3f') ' ( ' num2str((n+1)/N_t*100, '%4.1f') '% )']);

    % Plot q
%     f = figure(3);
%     f.Position = [100, 600, 1500, 400];
%     plot(x_axis, q_next);
%     xlim([0,len_x]);
%     ylim([-c0,c0]);

    % pause(0.1);

end